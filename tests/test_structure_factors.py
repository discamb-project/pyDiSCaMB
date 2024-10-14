import pydiscamb

import pytest


@pytest.mark.parametrize(
    "model",
    [
        pydiscamb.calculate_structure_factors_IAM,
        pydiscamb.calculate_structure_factors_TAAM,
    ],
)
class TestSF:
    def test_is_callable(self, model):
        assert callable(model)

    def test_basic_call(self, model, tyrosine):
        model(tyrosine, 3)

    def test_output(self, model, tyrosine):
        output = model(tyrosine, 5)
        assert hasattr(output, "__iter__")
        assert hasattr(output, "__getitem__")
        assert isinstance(output[0], complex)
