"""MoleCode — LLM system prompts.

The authoritative, human/LLM-readable specifications of the MoleCode Mermaid
grammar. Feed these as the system prompt when asking an LLM to read, write, or
edit molecules in MoleCode.

    MOLECULE_SYSTEM_PROMPT  — small-molecule grammar (atoms / bonds / stereo).
    MARKUSH_SYSTEM_PROMPT   — Markush extension (``{}`` abbreviation nodes).
"""

from .molecule_system_prompt import BASE_INSTRUCTION as MOLECULE_SYSTEM_PROMPT
from .markush_system_prompt import MARKUSH_SYSTEM_PROMPT

__all__ = ["MOLECULE_SYSTEM_PROMPT", "MARKUSH_SYSTEM_PROMPT"]
