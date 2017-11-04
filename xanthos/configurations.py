"""
Model configurations for Xanthos.

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

from xanthos.components import Components


def hargreaves_hejazi_simple(config):
    """
    Model configuration for the following:

    PET:                Hargreaves
    RUNOFF:             Hejazi
    ROUTING:            Simple


    :param config:      Configuration object generated from user-defined config.ini file
    :return:            Object containing all return values.
    """

    # instantiate hydro class
    c = Components(config)

    # load climate data
    c.load_climate()

    # load reference data
    c.load_reference()

    # PET calculation related
    c.setup_pet()

    # MSMC related data
    c.setup_msmc()

    # Flow related data
    c.setup_flow()

    # allocation
    c.setup_alloc()

    # spin up
    c.simulation(num_steps=config.SpinUp,
                 pet=True,
                 runoff=True,
                 runoff_step='month',
                 routing=True,
                 routing_step='month',
                 notify='Spin Up')

    # run model
    c.simulation(num_steps=config.nmonths,
                 pet=True,
                 runoff=True,
                 runoff_step='month',
                 routing=True,
                 routing_step='month',
                 notify='Simulation')

    # accessible water module
    c.accessible_water()

    # hydropower potential
    c.hydropower_potential()

    # hydropower actual
    c.hydropower_actual()

    # diagnostics
    c.diagnostics()

    # output simulation data
    c.output_simulation()

    # aggregate outputs
    c.aggregate_outputs()

    # create time series plots
    c.plots()

    return c


def hargreaves_abcd_simple(config):
    """
    Model configuration for the following:

    PET:                Hargreaves
    RUNOFF:             ABCD
    ROUTING:            Simple


    Spin-up is built in to the ABCD model and does not need to be calculated separately.

    :param config:      Configuration object generated from user-defined config.ini file
    :return:            Object containing all return values.
    """

    # instantiate hydro class
    c = Components(config)

    # load climate data
    c.load_climate()

    # load reference data
    c.load_reference()

    # PET calculation related
    c.setup_pet()

    # MSMC related data
    c.setup_msmc()

    # Flow related data
    c.setup_flow()

    # allocation
    c.setup_alloc()

    # run model
    c.simulation(num_steps=config.nmonths,
                 pet=True,
                 runoff=True,
                 runoff_step=None,
                 routing=True,
                 routing_step='month',
                 notify='Simulation')

    # accessible water module
    c.accessible_water()

    # hydropower potential
    c.hydropower_potential()

    # hydropower actual
    c.hydropower_actual()

    # diagnostics
    c.diagnostics()

    # output simulation data
    c.output_simulation()

    # aggregate outputs
    c.aggregate_outputs()

    # create time series plots
    c.plots()

    return c
