"""
Model configurations for Xanthos.

@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

from components import Components


def hargreaves_gwam_mrtm(config, spinup=True):
    """
    Model configuration for the following:

    PET:                Hargreaves
    RUNOFF:             Global Water Assessment Model (GWAM)
    ROUTING:            Modified River Transport Model (MRTM)


    :param config:      Configuration object generated from user-defined config.ini file
    :return:            Object containing all return values.
    """

    # instantiate hydro class
    c = Components(config)

    # spin up
    c.simulation(pet=True,
                 pet_num_steps=config.runoff_spinup,
                 pet_step='month',
                 runoff=True,
                 runoff_num_steps=config.runoff_spinup,
                 runoff_step='month',
                 routing=True,
                 routing_num_steps=config.routing_spinup,
                 routing_step='month',
                 notify='Spin Up')

    # run model
    c.simulation(pet=True,
                 pet_num_steps=config.nmonths,
                 pet_step='month',
                 runoff=True,
                 runoff_num_steps=config.nmonths,
                 runoff_step='month',
                 routing=True,
                 routing_num_steps=config.nmonths,
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


def hargreaves_abcd_mrtm(config):
    """
    Model configuration for the following:

    PET:                Hargreaves
    RUNOFF:             ABCD
    ROUTING:            Modified River Transport Model (MRTM)


    Spin-up is built in to the ABCD model and does not need to be calculated separately.

    :param config:      Configuration object generated from user-defined config.ini file
    :return:            Object containing all return values.
    """

    # instantiate hydro class
    c = Components(config)

    # run model
    c.simulation(pet=True,
                 pet_num_steps=config.nmonths,
                 pet_step='month',
                 runoff=True,
                 runoff_step=None,
                 routing=True,
                 routing_num_steps=config.nmonths,
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


def none_none_mrtm(config):
    """
    Model configuration for the following:

    PET:                None
    RUNOFF:             None
    ROUTING:            Modified River Transport Model (MRTM)


    Providing a custom runoff (Q) array to only run the routing module.

    :param config:      Configuration object generated from user-defined config.ini file
    :return:            Object containing all return values.
    """

    # instantiate hydro class
    c = Components(config)

    # spin up
    c.simulation(pet=False,
                 runoff=False,
                 routing=True,
                 routing_num_steps=config.routing_spinup,
                 routing_step='month',
                 notify='Spin Up')
    # run model
    c.simulation(pet=False,
                 runoff=False,
                 routing=True,
                 routing_num_steps=config.nmonths,
                 routing_step='month',
                 notify='Simulation')

    # output simulation data
    c.output_simulation()

    # aggregate outputs
    c.aggregate_outputs()


def none_abcd_mrtm(config):
    """
    Model configuration for the following:

    PET:                None - input provided by user
    RUNOFF:             ABCD
    ROUTING:            Modified River Transport Model (MRTM)


    Spin-up is built in to the ABCD model and does not need to be calculated separately.

    :param config:      Configuration object generated from user-defined config.ini file
    :return:            Object containing all return values.
    """

    # instantiate hydro class
    c = Components(config)

    # run model
    c.simulation(pet=False,
                 pet_step=None,
                 runoff=True,
                 runoff_step=None,
                 routing=False,
                 routing_num_steps=config.nmonths,
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


def pm_abcd_mrtm(config):
    """
    Model configuration for the following:

    PET:                Penman-Monteith
    RUNOFF:             ABCD
    ROUTING:            Modified River Transport Model (MRTM)


    Spin-up is built in to the ABCD model and does not need to be calculated separately.

    :param config:      Configuration object generated from user-defined config.ini file
    :return:            Object containing all return values.
    """

    # instantiate hydro class
    c = Components(config)

    # run model
    c.simulation(pet=True,
                 pet_step=None,
                 runoff=True,
                 runoff_step=None,
                 routing=True,
                 routing_num_steps=config.nmonths,
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