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
    assert client.model == "MiniMax-M3.1-Flash"
    assert PROVIDER_MODELS["minimax"]["models"] == (
        "MiniMax-M3.1-Flash",
        "MiniMax-M3",
    )


def test_minimax_model_metadata_is_current():
    metadata = PROVIDER_MODELS["minimax"]["metadata"]
    assert tuple(metadata) == ("MiniMax-M3.1-Flash", "MiniMax-M3")
    flash = metadata["MiniMax-M3.1-Flash"]
    assert flash["context_window"] == 1_000_000
    assert flash["max_output_tokens"] == 524_288
    assert flash["pricing_usd_per_million_tokens"] == {
        "input": 0.15, "output": 0.6, "cache_read": 0.015, "cache_write": None,
    }
    assert flash["thinking"] == ("adaptive",)
    assert flash["effort"] == ("low", "medium", "high", "xhigh", "max")
    assert flash["default_effort"] == "max"
    m3 = metadata["MiniMax-M3"]
    assert m3["pricing_usd_per_million_tokens"] == {
        "input": 0.3, "output": 1.2, "cache_read": 0.06, "cache_write": None,
    }
    assert m3["pricing_usd_per_million_tokens_above_512k"] == {
        "input": 0.6, "output": 2.4, "cache_read": 0.12, "cache_write": None,
    }
    assert m3["priority_price_multiplier"] == 1.5
    assert m3["thinking"] == ("adaptive", "disabled")
    assert m3["default_thinking_by_transport"] == {
        "openai": "adaptive",
        "anthropic": "disabled",
    }
    assert m3["effort_ignored"] is True
    assert m3["input_modalities"] == flash["input_modalities"] == (
        "text", "image", "video",
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
    assert client.model == "MiniMax-M3.1-Flash"


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
    assert captured["headers"]["X-api-key"] == "dummy"
    assert captured["payload"] == {
        "model": "MiniMax-M3.1-Flash",
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
        "thinking": {"type": "adaptive"},
        "system": "Return MoleCode.",
    }


def test_anthropic_messages_sends_positive_temperature():
    client = LLMClient(
        api_key="dummy",
        provider="minimax",
        transport="anthropic",
    )
    captured = {}
    response = {"content": [{"type": "text", "text": "ok"}]}
    with mock.patch.object(
        llm_mod.urllib.request,
        "urlopen",
        _capture_urlopen(captured, response),
    ):
        assert client.chat("hi", temperature=0.5) == "ok"

    assert captured["payload"]["temperature"] == 0.5


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


@pytest.mark.parametrize("transport", ["openai", "anthropic"])
@pytest.mark.parametrize("region,host", [
    ("global_en", "api.minimax.io"), ("cn_zh", "api.minimaxi.com"),
])
@pytest.mark.parametrize("source,expected", [
    ("default", "MiniMax-M3.1-Flash"),
    ("constructor", "MiniMax-M3"),
    ("environment", "MiniMax-M2.7"),
    ("request", "custom-model"),
])
def test_serialized_model_resolution(monkeypatch, transport, region, host,
                                     source, expected):
    kwargs = {"api_key": "dummy", "transport": transport,
              "provider": "minimax", "region": region}
    request_kwargs = {}
    if source == "constructor":
        monkeypatch.setenv("MOLECODE_MODEL", "saved-model")
        kwargs["model"] = expected
    elif source == "environment":
        monkeypatch.setenv("MOLECODE_MODEL", expected)
    elif source == "request":
        kwargs["model"] = "MiniMax-M2.7"
        request_kwargs["model"] = expected
    client = LLMClient(**kwargs)
    captured = {}
    response = ({"choices": [{"message": {"content": "ok"}}]}
                if transport == "openai"
                else {"content": [{"type": "text", "text": "ok"}]})
    with mock.patch.object(llm_mod.urllib.request, "urlopen",
                           _capture_urlopen(captured, response)):
        assert client.chat("hi", **request_kwargs) == "ok"
    suffix = "/v1/chat/completions" if transport == "openai" else "/anthropic/v1/messages"
    assert captured["url"] == "https://" + host + suffix
    assert captured["payload"]["model"] == expected
    for field in ("thinking", "reasoning_effort", "output_config",
                  "service_tier", "reasoning_split"):
        assert field not in captured["payload"]


@pytest.mark.parametrize("transport", ["openai", "anthropic"])
@pytest.mark.parametrize("effort", ["low", "medium", "high", "xhigh", "max"])
def test_flash_thinking_and_effort_reach_transport(transport, effort):
    client = LLMClient(api_key="dummy", provider="minimax", transport=transport)
    options = {"thinking": {"type": "adaptive"}}
    if transport == "openai":
        options["reasoning_effort"] = effort
    else:
        options["output_config"] = {"effort": effort}
    captured = {}
    response = ({"choices": [{"message": {"content": "ok"}}]}
                if transport == "openai"
                else {"content": [{"type": "text", "text": "ok"}]})
    with mock.patch.object(llm_mod.urllib.request, "urlopen",
                           _capture_urlopen(captured, response)):
        assert client.chat("hi", **options) == "ok"
    assert captured["payload"]["model"] == "MiniMax-M3.1-Flash"
    for field, value in options.items():
        assert captured["payload"][field] == value


@pytest.mark.parametrize("transport", ["openai", "anthropic"])
def test_saved_model_and_custom_endpoint_keep_request_options(transport):
    client = LLMClient(api_key="dummy", provider="minimax", transport=transport,
                       base_url="https://custom.example.com", model="MiniMax-M3")
    captured = {}
    response = ({"choices": [{"message": {"content": "ok"}}]}
                if transport == "openai"
                else {"content": [{"type": "text", "text": "ok"}]})
    with mock.patch.object(llm_mod.urllib.request, "urlopen",
                           _capture_urlopen(captured, response)):
        assert client.chat("hi", thinking={"type": "disabled"},
                           service_tier="priority") == "ok"
    assert captured["url"].startswith("https://custom.example.com/")
    assert captured["payload"]["model"] == "MiniMax-M3"
    assert captured["payload"]["thinking"] == {"type": "disabled"}
    assert captured["payload"]["service_tier"] == "priority"


@pytest.mark.parametrize("transport", ["openai", "anthropic"])
def test_call_llm_uses_environment_provider_default(monkeypatch, transport):
    monkeypatch.setenv("MOLECODE_PROVIDER", "minimax")
    monkeypatch.setenv("MOLECODE_TRANSPORT", transport)
    monkeypatch.setenv("MOLECODE_REGION", "cn_zh")
    captured = {}
    response = ({"choices": [{"message": {"content": "ok"}}]}
                if transport == "openai"
                else {"content": [{"type": "text", "text": "ok"}]})
    with mock.patch.object(llm_mod.urllib.request, "urlopen",
                           _capture_urlopen(captured, response)):
        assert llm_mod.call_llm("system", "hi", api_key="dummy") == "ok"
    assert captured["payload"]["model"] == "MiniMax-M3.1-Flash"
    assert captured["url"].startswith("https://api.minimaxi.com/")


def test_non_minimax_default_payload_is_unchanged():
    client = LLMClient(api_key="dummy")
    captured = {}
    response = {"choices": [{"message": {"content": "ok"}}]}
    with mock.patch.object(llm_mod.urllib.request, "urlopen",
                           _capture_urlopen(captured, response)):
        assert client.chat("hi") == "ok"
    assert captured["payload"] == {
        "model": llm_mod.DEFAULT_MODEL,
        "messages": [{"role": "user", "content": "hi"}],
        "temperature": 0.0,
    }


def test_invalid_transport_is_rejected():
    with pytest.raises(ValueError, match="Unsupported transport"):
        LLMClient(api_key="dummy", transport="unknown")


def test_anthropic_transport_requires_base_url_or_provider():
    with pytest.raises(ValueError, match="Anthropic transport requires"):
        LLMClient(api_key="dummy", transport="anthropic")
