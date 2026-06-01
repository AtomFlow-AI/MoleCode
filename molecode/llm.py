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
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from typing import Any, Dict, List, Optional

DEFAULT_BASE_URL = "https://api.openai.com/v1"
DEFAULT_MODEL = "gemini-3.1-pro-preview"

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
    model:
        Default model name. Falls back to ``$MOLECODE_MODEL`` then
        ``gemini-3.1-pro-preview``.
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
        self.base_url = (
            base_url or os.environ.get("MOLECODE_BASE_URL") or DEFAULT_BASE_URL
        ).rstrip("/")
        self.model = model or os.environ.get("MOLECODE_MODEL") or DEFAULT_MODEL
        self.timeout = timeout
        self.default_temperature = default_temperature

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
