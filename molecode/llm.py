"""A minimal, dependency-free LLM client for driving MoleCode tasks.

MoleCode itself never needs an LLM — it is a pure representation library. This
module is an *optional* convenience so you can run the understanding /
generation / editing / reasoning workflows without wiring up an SDK yourself.

``LLMClient`` speaks the **OpenAI Chat Completions** protocol over plain stdlib
``urllib`` (no third-party packages), so it works with any OpenAI-compatible
endpoint — OpenAI, Azure OpenAI, DeepSeek, Together, vLLM, Ollama, etc. You
supply the API key and base URL; nothing is hard-coded.

    from molecode.llm import LLMClient
    from molecode.prompts import MOLECULE_SYSTEM_PROMPT

    client = LLMClient(api_key="sk-...", base_url="https://api.openai.com/v1",
                       model="gemini-3.1-pro-preview")
    answer = client.chat("How many carbons are in this graph? ...",
                         system=MOLECULE_SYSTEM_PROMPT)

Credentials may also come from the environment (so you never commit a key):

    MOLECODE_API_KEY   (or OPENAI_API_KEY)   — required
    MOLECODE_BASE_URL  — default https://api.openai.com/v1
    MOLECODE_MODEL     — default gemini-3.1-pro-preview

Prefer the official ``openai`` SDK? You don't need this class at all — the
MoleCode prompts are plain strings, so pass them straight to
``openai.OpenAI().chat.completions.create(...)``.

Provider presets
----------------
``PROVIDER_PRESETS`` maps a provider name to its OpenAI Chat Completions base
URL per region, so callers can opt into a provider without hard-coding URLs::

    from molecode.llm import LLMClient, PROVIDER_PRESETS

    client = LLMClient(api_key="...", provider="minimax", region="global_en",
                       model="MiniMax-M3")
    client = LLMClient(api_key="...", provider="minimax", region="cn_zh",
                       model="MiniMax-M2.7")

A preset resolves ``base_url`` from ``provider`` + ``region`` (falling back to
``MOLECODE_BASE_URL``/``OPENAI_BASE_URL`` then the OpenAI default); an explicit
``base_url`` always wins. ``MOLECODE_PROVIDER``/``MOLECODE_REGION`` mirror the
constructor kwargs for environment-only configuration.
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from typing import Any, Dict, List, Optional

DEFAULT_BASE_URL = "https://api.openai.com/v1"
DEFAULT_MODEL = "gemini-3.1-pro-preview"

# Provider-specific OpenAI Chat Completions base URLs, keyed by provider then
# region. Each entry mirrors the provider's documented ``openai_base_url``.
PROVIDER_PRESETS: Dict[str, Dict[str, str]] = {
    "minimax": {
        "global_en": "https://api.minimax.io/v1",
        "cn_zh": "https://api.minimaxi.com/v1",
    },
}

_IMAGE_MIME = {
    ".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg",
    ".gif": "image/gif", ".webp": "image/webp", ".bmp": "image/bmp",
}


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
    """OpenAI-compatible chat client. You provide ``api_key`` and ``base_url``.

    Parameters
    ----------
    api_key:
        Bearer token. Falls back to ``$MOLECODE_API_KEY`` then ``$OPENAI_API_KEY``.
    base_url:
        Chat-completions base URL (without the ``/chat/completions`` suffix).
        Falls back to ``$MOLECODE_BASE_URL`` then ``https://api.openai.com/v1``.
        Ignored when ``provider`` resolves to a preset (an explicit ``base_url``
        always wins over a preset).
    model:
        Default model name. Falls back to ``$MOLECODE_MODEL`` then
        ``gemini-3.1-pro-preview``.
    provider:
        Optional provider name (e.g. ``"minimax"``). When set, the base URL is
        resolved from ``PROVIDER_PRESETS`` using ``region`` unless ``base_url``
        is given explicitly. Falls back to ``$MOLECODE_PROVIDER``.
    region:
        Region key for the provider preset (e.g. ``"global_en"`` or
        ``"cn_zh"``). Falls back to ``$MOLECODE_REGION`` then the preset's first
        region.
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
        timeout: float = 120.0,
        default_temperature: float = 0.0,
    ) -> None:
        self.api_key = (
            api_key
            or os.environ.get("MOLECODE_API_KEY")
            or os.environ.get("OPENAI_API_KEY")
        )
        if not self.api_key:
            raise ValueError(
                "No API key provided. Pass api_key=... or set the "
                "MOLECODE_API_KEY (or OPENAI_API_KEY) environment variable."
            )
        self.provider = provider or os.environ.get("MOLECODE_PROVIDER")
        self.region = region or os.environ.get("MOLECODE_REGION")
        self.base_url = self._resolve_base_url(base_url).rstrip("/")
        self.model = model or os.environ.get("MOLECODE_MODEL") or DEFAULT_MODEL
        self.timeout = timeout
        self.default_temperature = default_temperature

    def _resolve_base_url(self, base_url: Optional[str]) -> str:
        """Resolve the chat base URL from an explicit value, preset, or env."""
        if base_url:
            return base_url
        env_url = os.environ.get("MOLECODE_BASE_URL") or os.environ.get(
            "OPENAI_BASE_URL"
        )
        if env_url:
            return env_url
        if self.provider and self.provider.lower() in PROVIDER_PRESETS:
            regions = PROVIDER_PRESETS[self.provider.lower()]
            key = self.region if self.region and self.region in regions else None
            if key is None:
                key = next(iter(regions))
            return regions[key]
        return DEFAULT_BASE_URL

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
        messages: List[Dict[str, str]],
        *,
        model: Optional[str] = None,
        temperature: Optional[float] = None,
        **extra: Any,
    ) -> str:
        """Send a full ``messages`` list, return the assistant's text content."""
        payload: Dict[str, Any] = {
            "model": model or self.model,
            "messages": messages,
            "temperature": (
                self.default_temperature if temperature is None else temperature
            ),
            **extra,
        }
        request = urllib.request.Request(
            f"{self.base_url}/chat/completions",
            data=json.dumps(payload).encode("utf-8"),
            headers={
                "Authorization": f"Bearer {self.api_key}",
                "Content-Type": "application/json",
            },
            method="POST",
        )
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as resp:
                data = json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:  # surface the server's message
            detail = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"LLM API HTTP {exc.code}: {detail}") from exc
        except urllib.error.URLError as exc:
            raise RuntimeError(f"LLM API connection error: {exc.reason}") from exc

        return data["choices"][0]["message"]["content"]


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
    except ValueError:
        return None
    return client.chat(user, system=system, temperature=temperature)
