import taam_sf

import pytest


@pytest.mark.parametrize("model", [taam_sf.test_IAM, taam_sf.test_TAAM])
class TestSF:

    def test_is_callable(self, model):
        assert callable(model)

    def test_basic_call(self, model, lysosyme):
        model(lysosyme, 3)

    def test_output(self, model, lysosyme):
        output = model(lysosyme, 5)
        assert hasattr(output, "__iter__")
        assert hasattr(output, "__getitem__")
        assert isinstance(output[0], complex)
