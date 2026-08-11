"""Tests for the provider/region preset support in ``molecode.llm``.

These exercise URL/model resolution offline (no network) so they run anywhere.
"""

import json
import os
from unittest import mock

import pytest

from molecode import llm as llm_mod
from molecode.llm import LLMClient, PROVIDER_MODELS, PROVIDER_PRESETS


class _FakeResponse:
    def __init__(self, payload):
        self.payload = payload

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self):
        return json.dumps(self.payload).encode("utf-8")


def _capture_urlopen(captured, response):
    def fake_urlopen(request, timeout=None):
        captured["url"] = request.full_url
        captured["headers"] = request.headers
        captured["payload"] = json.loads(request.data.decode("utf-8"))
        return _FakeResponse(response)

    return fake_urlopen


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
        "MOLECODE_TRANSPORT",
        "ANTHROPIC_API_KEY",
        "ANTHROPIC_AUTH_TOKEN",
        "ANTHROPIC_BASE_URL",
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
    assert client.provider == "minimax"


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


def test_minimax_uses_current_default_model():
    client = LLMClient(api_key="dummy", provider="minimax")
    assert client.model == "MiniMax-M3"
    assert PROVIDER_MODELS["minimax"]["models"] == (
        "MiniMax-M3",
        "MiniMax-M2.7",
    )


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
    response = {"choices": [{"message": {"content": "ok"}}]}
    with mock.patch.object(
        llm_mod.urllib.request,
        "urlopen",
        _capture_urlopen(captured, response),
    ):
        out = client.chat("hi")
    assert out == "ok"
    assert captured["url"] == "https://api.minimax.io/v1/chat/completions"
    assert captured["headers"]["Authorization"] == "Bearer dummy"


def test_openai_image_request_path_is_unchanged():
    client = LLMClient(api_key="dummy", provider="minimax", model="MiniMax-M3")
    captured = {}
    response = {"choices": [{"message": {"content": "ok"}}]}
    with mock.patch.object(
        llm_mod.urllib.request,
        "urlopen",
        _capture_urlopen(captured, response),
    ):
        out = client.chat("describe", images=["https://x.example.com/a.png"])

    assert out == "ok"
    assert captured["payload"]["messages"][0]["content"] == [
        {"type": "text", "text": "describe"},
        {
            "type": "image_url",
            "image_url": {"url": "https://x.example.com/a.png"},
        },
    ]


def test_presets_contain_both_minimax_regions():
    assert "minimax" in PROVIDER_PRESETS
    regions = PROVIDER_PRESETS["minimax"]
    assert regions["global_en"] == {
        "openai": "https://api.minimax.io/v1",
        "anthropic": "https://api.minimax.io/anthropic",
    }
    assert regions["cn_zh"] == {
        "openai": "https://api.minimaxi.com/v1",
        "anthropic": "https://api.minimaxi.com/anthropic",
    }


@pytest.mark.parametrize(
    ("region", "base_url"),
    [
        ("global_en", "https://api.minimax.io/anthropic"),
        ("cn_zh", "https://api.minimaxi.com/anthropic"),
    ],
)
def test_anthropic_transport_uses_regional_endpoint(region, base_url):
    client = LLMClient(
        api_key="dummy",
        provider="minimax",
        region=region,
        transport="anthropic",
    )
    assert client.base_url == base_url
    assert client.model == "MiniMax-M3"


def test_anthropic_messages_request_and_response():
    client = LLMClient(
        api_key="dummy",
        provider="minimax",
        region="global_en",
        transport="anthropic",
    )
    captured = {}
    response = {
        "content": [
            {"type": "thinking", "thinking": "..."},
            {"type": "text", "text": "first"},
            {"type": "text", "text": " second"},
        ]
    }
    with mock.patch.object(
        llm_mod.urllib.request,
        "urlopen",
        _capture_urlopen(captured, response),
    ):
        out = client.chat(
            "describe",
            system="Return MoleCode.",
            images=["https://x.example.com/a.png"],
            max_tokens=512,
            thinking={"type": "adaptive"},
        )

    assert out == "first second"
    assert captured["url"] == "https://api.minimax.io/anthropic/v1/messages"
    assert captured["headers"]["Authorization"] == "Bearer dummy"
    assert captured["payload"] == {
        "model": "MiniMax-M3",
        "messages": [
            {
                "role": "user",
                "content": [
                    {"type": "text", "text": "describe"},
                    {
                        "type": "image",
                        "source": {
                            "type": "url",
                            "url": "https://x.example.com/a.png",
                        },
                    },
                ],
            }
        ],
        "max_tokens": 512,
        "temperature": 0.0,
        "thinking": {"type": "adaptive"},
        "system": "Return MoleCode.",
    }


def test_anthropic_converts_local_image_to_base64(tmp_path):
    image = tmp_path / "image.png"
    image.write_bytes(b"png")
    client = LLMClient(
        api_key="dummy",
        provider="minimax",
        transport="anthropic",
        model="MiniMax-M3",
    )
    captured = {}
    response = {"content": [{"type": "text", "text": "ok"}]}
    with mock.patch.object(
        llm_mod.urllib.request,
        "urlopen",
        _capture_urlopen(captured, response),
    ):
        assert client.chat("describe", images=[str(image)]) == "ok"

    source = captured["payload"]["messages"][0]["content"][1]["source"]
    assert source == {"type": "base64", "media_type": "image/png", "data": "cG5n"}


def test_anthropic_base_url_environment_wins_over_preset():
    with mock.patch.dict(
        os.environ,
        {
            "MOLECODE_PROVIDER": "minimax",
            "MOLECODE_TRANSPORT": "anthropic",
            "ANTHROPIC_BASE_URL": "https://anthropic.example.com",
        },
    ):
        client = LLMClient(api_key="dummy")
    assert client.base_url == "https://anthropic.example.com"


def test_invalid_transport_is_rejected():
    with pytest.raises(ValueError, match="Unsupported transport"):
        LLMClient(api_key="dummy", transport="unknown")


def test_anthropic_transport_requires_base_url_or_provider():
    with pytest.raises(ValueError, match="Anthropic transport requires"):
        LLMClient(api_key="dummy", transport="anthropic")
