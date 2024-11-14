import pytest
import numpy as np

from cctbx.array_family import flex
import mmtbx.f_model

from pydiscamb import DiscambWrapper

def test_target_gradients(random_structure):

    # Use f_calc as f_obs
    d_min = 2
    f_obs = abs(random_structure.structure_factors(d_min=d_min).f_calc())

    # Set up model
    random_structure.shake_sites_in_place(rms_difference=0.1)
    scatterers = random_structure.scatterers()
    model = mmtbx.f_model.manager(f_obs=f_obs, xray_structure=random_structure)
    target = model.target_functor()(compute_gradients=True)

    # Get target derivatives
    d_target_d_fcalc = target.d_target_d_f_calc_work()

    # Calculate derivatives with discamb
    w = DiscambWrapper(random_structure)
    w.set_indices(d_target_d_fcalc.indices())
    discamb_result = w.d_target_d_params(list(d_target_d_fcalc.data()))

    ### Check equality
    iselection = flex.bool(scatterers.size(), True).iselection()

    # Site
    scatterers.flags_set_grads(state=False)
    scatterers.flags_set_grad_site(iselection=iselection)
    g = flex.vec3_double(target.gradients_wrt_atomic_parameters().packed())
    expected = np.array(g)
    actual = np.array([res.site_derivatives for res in discamb_result])
    assert pytest.approx(expected) == actual

    # ADP
    scatterers.flags_set_grads(state=False)
    scatterers.flags_set_grad_u_iso(iselection=iselection)
    g = target.gradients_wrt_atomic_parameters().packed()
    assert pytest.approx(np.array(g)) == np.array([res.adp_derivatives for res in discamb_result]).flatten()

    # Occupancy
    scatterers.flags_set_grads(state=False)
    scatterers.flags_set_grad_occupancy(iselection=iselection)
    g = target.gradients_wrt_atomic_parameters().packed()
    assert pytest.approx(list(g)) == [res.occupancy_derivatives for res in discamb_result]

def test_target_gradients_u_aniso(random_structure_u_aniso):

    # Use f_calc as f_obs
    d_min = 2
    f_obs = abs(random_structure_u_aniso.structure_factors(d_min=d_min).f_calc())

    # Set up model
    random_structure_u_aniso.shake_sites_in_place(rms_difference=0.1)
    scatterers = random_structure_u_aniso.scatterers()
    model = mmtbx.f_model.manager(f_obs=f_obs, xray_structure=random_structure_u_aniso)
    target = model.target_functor()(compute_gradients=True)

    # Get target derivatives
    d_target_d_fcalc = target.d_target_d_f_calc_work()

    # Calculate derivatives with discamb
    w = DiscambWrapper(random_structure_u_aniso)
    w.set_indices(d_target_d_fcalc.indices())
    discamb_result = w.d_target_d_params(list(d_target_d_fcalc.data()))

    ### Check equality
    iselection = flex.bool(scatterers.size(), True).iselection()
    scatterers.flags_set_grads(state=False)
    scatterers.flags_set_grad_u_aniso(iselection=iselection)
    g = target.gradients_wrt_atomic_parameters().packed()
    assert pytest.approx(np.array(g)) == np.array([res.adp_derivatives for res in discamb_result]).flatten()
