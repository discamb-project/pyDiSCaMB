from typing import TYPE_CHECKING, Dict

from pydiscamb.discamb_wrapper.discamb_wrapper import (
    DiscambWrapper,
    _concat_scatterer_labels,
)
from pydiscamb.discamb_wrapper.fcalc_method import FCalcMethod

if TYPE_CHECKING:
    from cctbx.xray.structure import structure


class DiscambWrapperCached(DiscambWrapper):

    __cache = {}

    def __init__(
        self, xrs: "structure", method: FCalcMethod = FCalcMethod.IAM, **kwargs
    ):
        if self.__check_cache(xrs, method, kwargs) is not None:
            self.update_structure(xrs)
            return

        super().__init__(xrs, method, **kwargs)

        self.__cache[self.__get_cache_key(xrs, method, kwargs)] = self

    def __new__(cls, xrs: "structure", method: FCalcMethod = FCalcMethod.IAM, **kwargs):
        cache = cls.__check_cache(xrs, method, kwargs)
        if cache is None:
            return super().__new__(cls)
        return cache

    @classmethod
    def __check_cache(
        cls, xrs: "structure", method: FCalcMethod, kwargs: dict[str, str]
    ):
        key = cls.__get_cache_key(xrs, method, kwargs)
        return cls.__cache.get(key)

    @classmethod
    def __get_cache_key(
        cls, xrs: "structure", method: FCalcMethod, kwargs: Dict[str, str]
    ):
        atomstr = _concat_scatterer_labels(xrs)
        unitcell = xrs.unit_cell()
        params = method.to_cache_lookup_key(xrs, kwargs)
        key = (atomstr, unitcell, *params)
        print(key)
        return key
