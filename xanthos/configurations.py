"""
Model configurations for Xanthos.

@author:  Chris R. Vernon, Caleb J. Braun
@email:   chris.vernon@pnnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2019, Battelle Memorial Institute
"""

from xanthos.components import Components
import logging


class ConfigRunner:
    """
    Run the components specified by the configuration file.

    Model configurations are implemented for the for the following:

    PET:         One of:
                        'hs'            =  Hargreaves-Samani
                        'hargreaves'    =  Hargreaves
                        'pm'            =  Penman-Monteith
                        'thornthwaite'  =  Thornthwaite
    RUNOFF:      One of:
                        'abcd'  =  The ABCD model
                        'gwam'  =  Global Water Assessment Model
    ROUTING:     One of:
                        'mrtm'  =  Modified River Transport Model
                        'mrtm_managed' = MRTM with water mangement

    Available output statistics modules include:
        - Hydropower actual
        - Hydropower potential
        - Drought
        - Diagnostics

    Spin-up is built in to the ABCD model and does not need to be calculated
    separately. If any component needs to have the model spin up before
    running, they must be in the SPINUP_COMPONENTS list.

    It is possible to run any subset of the three main modules, as long as all
    the required inputs are provided. For example, to run only the routing
    module, provide a custom runoff array.

    Xanthos suports two types of modules: those that iterate one timestep at a
    time (ex. GWAM), and those that take all times at once (ex. ABCD). These
    properties must be specified in the *_STEPWISE_COMPONENTS in order to work
    properly.
    """

    def __init__(self, config):
        """
        Verify that the given configuration contains valid options.

        :param config:      Configuration object generated from user-defined config.ini file
        """
        # acceptable components to run
        PET_COMPONENTS = ['hs', 'hargreaves', 'pm', 'thornthwaite']
        RUNOFF_COMPONENTS = ['abcd', 'gwam', 'abcd_managed']
        ROUTING_COMPONENTS = ['mrtm', 'mrtm_managed']

        # runoff components that need the whole model to do spin-up
        SPINUP_COMPONENTS = ['gwam']

        # components that need to be run step-by-step (iterates over columns)
        PET_STEPWISE_COMPONENTS = ['hargreaves']
        RUNOFF_STEPWISE_COMPONENTS = ['gwam']
        ROUTING_STEPWISE_COMPONENTS = []

        # only run verified components
        self.run_pet = config.pet_module in PET_COMPONENTS
        self.run_runoff = config.runoff_module in RUNOFF_COMPONENTS
        self.run_routing = config.routing_module in ROUTING_COMPONENTS

        # check if we need to run spinup
        self.spinup = config.runoff_module in SPINUP_COMPONENTS

        # determine whether to run each module stepwise (iterate over each month), or
        # if they handle it internally (indicated by passing a value of 0)
        self.pet_timestep = config.nmonths * (config.pet_module in PET_STEPWISE_COMPONENTS)
        self.runoff_timestep = config.nmonths * (config.runoff_module in RUNOFF_STEPWISE_COMPONENTS)
        self.routing_timestep = config.nmonths * (config.routing_module in ROUTING_STEPWISE_COMPONENTS)

        self.config = config

    def run(self):
        """
        Run all valid components.

        :return:   Components object containing all return values.
        """
        # make sure there is something to run
        if not (self.run_pet or self.run_runoff or self.run_routing):
            logging.warning("Selected configuration {0} not supported.".format(self.config.mod_cfg))
            return

        # instantiate variables for all components
        c = Components(self.config)

        if self.spinup:
            c.simulation(run_pet=self.run_pet,
                         run_runoff=self.run_runoff,
                         run_routing=self.run_routing,
                         pet_num_steps=self.config.runoff_spinup,
                         runoff_num_steps=self.config.runoff_spinup,
                         routing_num_steps=self.config.routing_spinup,
                         notify='Spin Up')

        c.simulation(run_pet=self.run_pet,
                     run_runoff=self.run_runoff,
                     run_routing=self.run_routing,
                     pet_num_steps=self.pet_timestep,
                     runoff_num_steps=self.runoff_timestep,
                     routing_num_steps=self.routing_timestep,
                     notify='Simulation')

        # accessible water module
        c.accessible_water()

        # drought module
        c.drought()

        # hydropower potential
        c.hydropower_potential()

        # hydropower actual
        c.hydropower_actual()

        # diagnostics
        c.diagnostics()

        # output simulation data
        c.output_simulation()

        # create time series plots
        c.plots()

        return c
