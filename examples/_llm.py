"""Thin example helper around the MoleCode LLM client.

The real client lives in the library: ``molecode.llm.LLMClient`` (an
OpenAI-compatible, dependency-free chat client). This module only adds the
"dry run" behaviour the examples want: if no API key is configured, print the
exact prompt that *would* be sent and return ``None`` so the script still runs
offline.

To actually call a model, configure credentials via environment variables
(any OpenAI-compatible endpoint works):

    MOLECODE_API_KEY    your key            (required to actually call)
    MOLECODE_BASE_URL   chat base URL       (default https://api.openai.com/v1)
    MOLECODE_MODEL      model name          (default gpt-4o-mini)

…or construct the client directly in your own code:

    from molecode.llm import LLMClient
    client = LLMClient(api_key="sk-...", base_url="...", model="...")
    reply = client.chat(user_prompt, system=system_prompt)
"""

from __future__ import annotations

from molecode.llm import call_llm as _call_llm


def call_llm(system: str, user: str, *, temperature: float = 0.0) -> str | None:
    """Call the LLM, or dry-run (print the prompt and return None) if no key."""
    reply = _call_llm(system, user, temperature=temperature)
    if reply is None:
        print("=" * 70)
        print("DRY RUN — no MOLECODE_API_KEY/OPENAI_API_KEY set.")
        print("Below is the exact prompt MoleCode would send to any LLM:\n")
        print("----- SYSTEM -----\n" + system)
        print("\n----- USER -----\n" + user)
        print("=" * 70)
    return reply
