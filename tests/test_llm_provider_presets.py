"""Tests for the provider/region preset support in ``molecode.llm``.

These exercise URL/model resolution offline (no network) so they run anywhere.
"""

import json
import os
from unittest import mock

import pytest

from molecode import llm as llm_mod
from molecode.llm import LLMClient, PROVIDER_PRESETS


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    """Ensure provider/base_url/model env vars do not leak between tests."""
    for name in (
        "MOLECODE_API_KEY",
        "OPENAI_API_KEY",
        "MOLECODE_BASE_URL",
        "OPENAI_BASE_URL",
        "MOLECODE_MODEL",
        "MOLECODE_PROVIDER",
        "MOLECODE_REGION",
    ):
        monkeypatch.delenv(name, raising=False)
    yield


def test_minimax_preset_global_endpoint():
    client = LLMClient(
        api_key="dummy",
        provider="minimax",
        region="global_en",
        model="MiniMax-M3",
    )
    assert client.base_url == "https://api.minimax.io/v1"
    assert client.model == "MiniMax-M3"
    assert client.provider == "minimax"
    assert client.region == "global_en"


def test_minimax_preset_china_endpoint():
    client = LLMClient(
        api_key="dummy",
        provider="minimax",
        region="cn_zh",
        model="MiniMax-M2.7",
    )
    assert client.base_url == "https://api.minimaxi.com/v1"
    assert client.model == "MiniMax-M2.7"


def test_provider_is_case_insensitive():
    client = LLMClient(api_key="dummy", provider="MiniMax", model="MiniMax-M3")
    assert client.base_url == "https://api.minimax.io/v1"


def test_explicit_base_url_wins_over_preset():
    client = LLMClient(
        api_key="dummy",
        provider="minimax",
        base_url="https://custom.example.com/v1",
        model="MiniMax-M3",
    )
    assert client.base_url == "https://custom.example.com/v1"


def test_unknown_region_falls_back_to_first():
    client = LLMClient(
        api_key="dummy",
        provider="minimax",
        region="mars",
        model="MiniMax-M3",
    )
    # First region in the preset dict is global_en.
    assert client.base_url == "https://api.minimax.io/v1"


def test_no_provider_uses_openai_default():
    client = LLMClient(api_key="dummy")
    assert client.base_url == "https://api.openai.com/v1"


def test_provider_via_environment():
    with mock.patch.dict(
        os.environ,
        {"MOLECODE_PROVIDER": "minimax", "MOLECODE_REGION": "cn_zh"},
    ):
        client = LLMClient(api_key="dummy", model="MiniMax-M3")
        assert client.base_url == "https://api.minimaxi.com/v1"


def test_base_url_env_wins_over_provider_env():
    with mock.patch.dict(
        os.environ,
        {
            "MOLECODE_PROVIDER": "minimax",
            "MOLECODE_BASE_URL": "https://env.example.com/v1",
        },
    ):
        client = LLMClient(api_key="dummy", model="MiniMax-M3")
        assert client.base_url == "https://env.example.com/v1"


def test_chat_completions_url_assembles_from_preset():
    client = LLMClient(
        api_key="dummy",
        provider="minimax",
        region="global_en",
        model="MiniMax-M3",
    )
    # complete() posts to "<base_url>/chat/completions"; verify without a
    # network call by stubbing urlopen.
    captured = {}

    class _FakeResp:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def read(self):
            return json.dumps({"choices": [{"message": {"content": "ok"}}]}).encode("utf-8")

    def _fake_urlopen(req, timeout=None):
        captured["url"] = req.full_url
        captured["headers"] = req.headers
        return _FakeResp()

    with mock.patch.object(llm_mod.urllib.request, "urlopen", _fake_urlopen):
        out = client.chat("hi")
    assert out == "ok"
    assert captured["url"] == "https://api.minimax.io/v1/chat/completions"
    assert captured["headers"]["Authorization"] == "Bearer dummy"


def test_presets_contain_both_minimax_regions():
    assert "minimax" in PROVIDER_PRESETS
    regions = PROVIDER_PRESETS["minimax"]
    assert regions["global_en"] == "https://api.minimax.io/v1"
    assert regions["cn_zh"] == "https://api.minimaxi.com/v1"


def test_image_path_still_supported_with_provider():
    # The image-capable request path must keep working alongside presets.
    client = LLMClient(
        api_key="dummy",
        provider="minimax",
        region="global_en",
        model="MiniMax-M3",
    )
    assert client.base_url == "https://api.minimax.io/v1"
    # image_to_data_uri is the helper used by chat(images=...).
    assert (
        llm_mod.image_to_data_uri("https://x.example.com/a.png")
        == "https://x.example.com/a.png"
    )
