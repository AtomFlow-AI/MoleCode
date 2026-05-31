"""Tiny optional LLM helper for the task examples.

The task examples (understanding / generation / editing / reasoning) build the
prompt that you would send to *any* LLM. To stay dependency-free and runnable
without an API key, ``call_llm`` works in two modes:

* If ``MOLECODE_API_KEY`` (or ``OPENAI_API_KEY``) is set, it sends the prompt to
  an OpenAI-compatible chat endpoint and returns the model's reply.
  Configure with environment variables:

      MOLECODE_API_KEY    your key            (required to actually call)
      MOLECODE_BASE_URL   chat completions base URL
                          (default https://api.openai.com/v1)
      MOLECODE_MODEL      model name          (default gpt-4o-mini)

* Otherwise it prints the assembled system+user prompt and returns ``None`` so
  the example still runs end-to-end (offline "dry run").

MoleCode itself never calls an LLM — this helper exists only so the examples are
copy-pasteable into your own pipeline.
"""

from __future__ import annotations

import json
import os
import urllib.request


def call_llm(system: str, user: str, *, temperature: float = 0.0) -> str | None:
    """Send ``system``/``user`` to an OpenAI-compatible endpoint, or dry-run."""
    api_key = os.environ.get("MOLECODE_API_KEY") or os.environ.get("OPENAI_API_KEY")
    if not api_key:
        print("=" * 70)
        print("DRY RUN — no MOLECODE_API_KEY/OPENAI_API_KEY set.")
        print("Below is the exact prompt you would send to any LLM:\n")
        print("----- SYSTEM -----\n" + system)
        print("\n----- USER -----\n" + user)
        print("=" * 70)
        return None

    base = os.environ.get("MOLECODE_BASE_URL", "https://api.openai.com/v1").rstrip("/")
    model = os.environ.get("MOLECODE_MODEL", "gpt-4o-mini")
    payload = {
        "model": model,
        "temperature": temperature,
        "messages": [
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
    }
    req = urllib.request.Request(
        f"{base}/chat/completions",
        data=json.dumps(payload).encode(),
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
    )
    with urllib.request.urlopen(req) as resp:
        data = json.loads(resp.read())
    return data["choices"][0]["message"]["content"]
