"""A minimal, dependency-free LLM client for driving MoleCode tasks.

MoleCode itself never needs an LLM — it is a pure representation library. This
module is an *optional* convenience so you can run the understanding /
generation / editing / reasoning workflows without wiring up an SDK yourself.

``LLMClient`` speaks the **OpenAI Chat Completions** and **Anthropic Messages**
protocols over plain stdlib ``urllib`` (no third-party packages). You supply the
API key and can either choose a provider preset or pass a base URL directly.

    from molecode.llm import LLMClient
    from molecode.prompts import MOLECULE_SYSTEM_PROMPT

    client = LLMClient(api_key="sk-...", base_url="https://api.openai.com/v1",
                       model="gemini-3.1-pro-preview")
    answer = client.chat("How many carbons are in this graph? ...",
                         system=MOLECULE_SYSTEM_PROMPT)

Credentials may also come from the environment (so you never commit a key):

    MOLECODE_API_KEY   (or OPENAI_API_KEY)   — required
    MOLECODE_BASE_URL  - explicit transport base URL
    MOLECODE_MODEL     - model override
    MOLECODE_PROVIDER  - provider preset name
    MOLECODE_REGION    - provider region
    MOLECODE_TRANSPORT - openai (default) or anthropic

Prefer the official ``openai`` SDK? You don't need this class at all — the
MoleCode prompts are plain strings, so pass them straight to
``openai.OpenAI().chat.completions.create(...)``.

Provider presets
----------------
``PROVIDER_PRESETS`` maps a provider and region to OpenAI Chat Completions and
Anthropic Messages base URLs, so callers can opt into either transport without
hard-coding URLs::

    from molecode.llm import LLMClient, PROVIDER_PRESETS

    client = LLMClient(api_key="...", provider="minimax", region="global_en",
                       model="MiniMax-M3")
    client = LLMClient(api_key="...", provider="minimax", region="cn_zh",
                       transport="anthropic", model="MiniMax-M2.7")

A preset resolves ``base_url`` from ``provider`` + ``region`` + ``transport``.
An explicit ``base_url`` always wins, followed by the matching environment
variable. The provider also supplies its current default model when ``model``
and ``MOLECODE_MODEL`` are unset.
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from typing import Any, Dict, List, Optional

DEFAULT_BASE_URL = "https://api.openai.com/v1"
DEFAULT_MODEL = "gemini-3.1-pro-preview"
SUPPORTED_TRANSPORTS = ("openai", "anthropic")

# Provider-specific API base URLs, keyed by provider, region, then transport.
PROVIDER_PRESETS: Dict[str, Dict[str, Dict[str, str]]] = {
    "minimax": {
        "global_en": {
            "openai": "https://api.minimax.io/v1",
            "anthropic": "https://api.minimax.io/anthropic",
        },
        "cn_zh": {
            "openai": "https://api.minimaxi.com/v1",
            "anthropic": "https://api.minimaxi.com/anthropic",
        },
    },
}

PROVIDER_MODELS: Dict[str, Dict[str, Any]] = {
    "minimax": {
        "default": "MiniMax-M3",
        "models": ("MiniMax-M3", "MiniMax-M2.7"),
        "metadata": {
            "MiniMax-M3": {
                "context_window": 1_000_000,
                "pricing_usd_per_million_tokens": {
                    "input": 0.6,
                    "output": 2.4,
                    "cache_read": 0.12,
                    "cache_write": None,
                },
                "input_modalities": ("text", "image", "video"),
                "thinking": ("adaptive", "disabled"),
            },
            "MiniMax-M2.7": {
                "context_window": 204_800,
                "pricing_usd_per_million_tokens": {
                    "input": 0.3,
                    "output": 1.2,
                    "cache_read": 0.06,
                    "cache_write": 0.375,
                },
                "input_modalities": ("text",),
                "thinking": ("always_on",),
            },
        },
    },
}

_IMAGE_MIME = {
    ".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg",
    ".gif": "image/gif", ".webp": "image/webp", ".bmp": "image/bmp",
}


class MissingAPIKeyError(ValueError):
    """Raised when no API key is available for an LLM request."""


def image_to_data_uri(path_or_url: str) -> str:
    """Return a value suitable for an OpenAI ``image_url`` content block.

    A remote ``http(s)://`` URL is returned unchanged; a local file path is read
    and encoded as a base64 ``data:`` URI (mime guessed from the extension).
    """
    import base64
    import os.path

    if path_or_url.startswith(("http://", "https://", "data:")):
        return path_or_url
    ext = os.path.splitext(path_or_url)[1].lower()
    mime = _IMAGE_MIME.get(ext, "image/png")
    with open(path_or_url, "rb") as fh:
        b64 = base64.b64encode(fh.read()).decode("ascii")
    return f"data:{mime};base64,{b64}"


class LLMClient:
    """Minimal OpenAI Chat Completions and Anthropic Messages client.

    Parameters
    ----------
    api_key:
        Bearer token. Falls back to ``$MOLECODE_API_KEY`` then ``$OPENAI_API_KEY``.
    base_url:
        Transport base URL without ``/chat/completions`` or ``/v1/messages``.
        An explicit value wins over environment variables and provider presets.
    model:
        Default model name. Falls back to ``$MOLECODE_MODEL``, the provider's
        default model, then ``gemini-3.1-pro-preview``.
    provider:
        Optional provider name (e.g. ``"minimax"``). When set, the base URL is
        resolved from ``PROVIDER_PRESETS`` using ``region`` unless ``base_url``
        is given explicitly. Falls back to ``$MOLECODE_PROVIDER``.
    region:
        Region key for the provider preset (e.g. ``"global_en"`` or
        ``"cn_zh"``). Falls back to ``$MOLECODE_REGION`` then the preset's first
        region.
    transport:
        ``"openai"`` (default) or ``"anthropic"``. Falls back to
        ``$MOLECODE_TRANSPORT``.
    timeout:
        Per-request timeout in seconds.
    default_temperature:
        Temperature used when ``chat``/``complete`` don't override it.
    """

    def __init__(
        self,
        api_key: Optional[str] = None,
        base_url: Optional[str] = None,
        model: Optional[str] = None,
        *,
        provider: Optional[str] = None,
        region: Optional[str] = None,
        transport: Optional[str] = None,
        timeout: float = 120.0,
        default_temperature: float = 0.0,
    ) -> None:
        self.transport = (
            transport or os.environ.get("MOLECODE_TRANSPORT") or "openai"
        ).lower()
        if self.transport not in SUPPORTED_TRANSPORTS:
            supported = ", ".join(SUPPORTED_TRANSPORTS)
            raise ValueError(
                f"Unsupported transport {self.transport!r}; expected {supported}."
            )
        self.api_key = (
            api_key
            or os.environ.get("MOLECODE_API_KEY")
            or (
                os.environ.get("ANTHROPIC_API_KEY")
                if self.transport == "anthropic"
                else None
            )
            or (
                os.environ.get("ANTHROPIC_AUTH_TOKEN")
                if self.transport == "anthropic"
                else None
            )
            or os.environ.get("OPENAI_API_KEY")
        )
        if not self.api_key:
            raise MissingAPIKeyError(
                "No API key provided. Pass api_key=... or set the "
                "MOLECODE_API_KEY, OPENAI_API_KEY, ANTHROPIC_API_KEY, or "
                "ANTHROPIC_AUTH_TOKEN environment variable."
            )
        configured_provider = provider or os.environ.get("MOLECODE_PROVIDER")
        self.provider = configured_provider.lower() if configured_provider else None
        self.region = region or os.environ.get("MOLECODE_REGION")
        self.base_url = self._resolve_base_url(base_url).rstrip("/")
        self.model = (
            model
            or os.environ.get("MOLECODE_MODEL")
            or self._provider_default_model()
            or DEFAULT_MODEL
        )
        self.timeout = timeout
        self.default_temperature = default_temperature

    def _resolve_base_url(self, base_url: Optional[str]) -> str:
        """Resolve the chat base URL from an explicit value, preset, or env."""
        if base_url:
            return base_url
        transport_env = (
            "ANTHROPIC_BASE_URL"
            if self.transport == "anthropic"
            else "OPENAI_BASE_URL"
        )
        env_url = os.environ.get("MOLECODE_BASE_URL") or os.environ.get(transport_env)
        if env_url:
            return env_url
        if self.provider and self.provider in PROVIDER_PRESETS:
            regions = PROVIDER_PRESETS[self.provider]
            key = self.region if self.region and self.region in regions else None
            if key is None:
                key = next(iter(regions))
            return regions[key][self.transport]
        if self.transport == "anthropic":
            raise ValueError(
                "Anthropic transport requires base_url, ANTHROPIC_BASE_URL, "
                "or a provider preset."
            )
        return DEFAULT_BASE_URL

    def _provider_default_model(self) -> Optional[str]:
        if not self.provider:
            return None
        provider_models = PROVIDER_MODELS.get(self.provider)
        if not provider_models:
            return None
        return provider_models["default"]

    def chat(
        self,
        user: str,
        system: Optional[str] = None,
        *,
        images: Optional[List[str]] = None,
        model: Optional[str] = None,
        temperature: Optional[float] = None,
        **extra: Any,
    ) -> str:
        """Single-turn helper: send a ``system`` + ``user`` message, return text.

        Pass ``images`` (a list of local file paths or URLs) to send a
        multimodal request to a vision-capable model — used for OCSR
        (molecule image -> MoleCode graph). Requires a model that accepts image
        input (e.g. gpt-4o, gpt-4o-mini, gemini, claude vision models).
        """
        messages: List[Dict[str, Any]] = []
        if system:
            messages.append({"role": "system", "content": system})
        if images:
            content: List[Dict[str, Any]] = [{"type": "text", "text": user}]
            for img in images:
                content.append({
                    "type": "image_url",
                    "image_url": {"url": image_to_data_uri(img)},
                })
            messages.append({"role": "user", "content": content})
        else:
            messages.append({"role": "user", "content": user})
        return self.complete(messages, model=model, temperature=temperature, **extra)

    def complete(
        self,
        messages: List[Dict[str, Any]],
        *,
        model: Optional[str] = None,
        temperature: Optional[float] = None,
        **extra: Any,
    ) -> str:
        """Send a full ``messages`` list, return the assistant's text content."""
        if self.transport == "anthropic":
            return self._complete_anthropic(
                messages,
                model=model,
                temperature=temperature,
                **extra,
            )

        payload: Dict[str, Any] = {
            "model": model or self.model,
            "messages": messages,
            "temperature": (
                self.default_temperature if temperature is None else temperature
            ),
            **extra,
        }
        data = self._post_json(
            f"{self.base_url}/chat/completions",
            payload,
            {
                "Authorization": f"Bearer {self.api_key}",
                "Content-Type": "application/json",
            },
        )
        return data["choices"][0]["message"]["content"]

    def _complete_anthropic(
        self,
        messages: List[Dict[str, Any]],
        *,
        model: Optional[str],
        temperature: Optional[float],
        **extra: Any,
    ) -> str:
        system_parts: List[str] = []
        converted_messages: List[Dict[str, Any]] = []
        for message in messages:
            role = message["role"]
            content = message["content"]
            if role == "system":
                if isinstance(content, str):
                    system_parts.append(content)
                else:
                    system_parts.extend(
                        block["text"]
                        for block in content
                        if block.get("type") == "text"
                    )
                continue
            converted_messages.append(
                {"role": role, "content": self._anthropic_content(content)}
            )

        payload: Dict[str, Any] = {
            "model": model or self.model,
            "messages": converted_messages,
            "temperature": (
                self.default_temperature if temperature is None else temperature
            ),
            **extra,
        }
        if system_parts:
            payload["system"] = "\n\n".join(system_parts)

        data = self._post_json(
            f"{self.base_url}/v1/messages",
            payload,
            {
                "Authorization": f"Bearer {self.api_key}",
                "Content-Type": "application/json",
            },
        )
        return "".join(
            block["text"]
            for block in data["content"]
            if block.get("type") == "text"
        )

    @staticmethod
    def _anthropic_content(content: Any) -> Any:
        if isinstance(content, str):
            return content

        converted: List[Dict[str, Any]] = []
        for block in content:
            if block.get("type") != "image_url":
                converted.append(block)
                continue
            image_url = block["image_url"]
            value = image_url["url"] if isinstance(image_url, dict) else image_url
            if value.startswith("data:"):
                metadata, data = value.split(",", 1)
                media_type = metadata[5:].split(";", 1)[0]
                source = {
                    "type": "base64",
                    "media_type": media_type,
                    "data": data,
                }
            else:
                source = {"type": "url", "url": value}
            converted.append({"type": "image", "source": source})
        return converted

    def _post_json(
        self,
        url: str,
        payload: Dict[str, Any],
        headers: Dict[str, str],
    ) -> Dict[str, Any]:
        request = urllib.request.Request(
            url,
            data=json.dumps(payload).encode("utf-8"),
            headers=headers,
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:  # surface the server's message
            detail = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"LLM API HTTP {exc.code}: {detail}") from exc
        except urllib.error.URLError as exc:
            raise RuntimeError(f"LLM API connection error: {exc.reason}") from exc


def call_llm(
    system: str,
    user: str,
    *,
    temperature: float = 0.0,
    **client_kwargs: Any,
) -> Optional[str]:
    """Convenience wrapper used by the examples.

    Constructs an :class:`LLMClient` from arguments/environment and returns the
    model reply. If no API key is configured, returns ``None`` instead of
    raising, so example scripts can "dry run" and just print the prompt.
    """
    try:
        client = LLMClient(**client_kwargs)
    except MissingAPIKeyError:
        return None
    return client.chat(user, system=system, temperature=temperature)
