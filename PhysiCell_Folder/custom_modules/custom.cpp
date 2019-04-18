/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

// declare cell definitions here 

Cell_Definition wound_cell; 
Cell_Definition bacterial_colony;
Cell_Definition bacterial_colony2;

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "default cell"; 
	
	// Setting Cell Cycle 
	cell_defaults.functions.cycle_model = live; 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	//int ncycle_Start_i = live.find_phase_index( PhysiCell_constants::live );
	//int ncycle_End_i = live.find_phase_index( PhysiCell_constants::live );
	int ncycle_Start_i = live.find_phase_index( PhysiCell_constants::live );
	

		// No Proliferation, Apoptosis, and Necrosis
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 
	//cell_defaults.phenotype.cycle.data.transition_rate(ncycle_Start_i,ncycle_End_i) = 0.0;
	
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0;


	// No phenotype change; 
	cell_defaults.functions.update_phenotype = NULL; 

	
	// Setting Orientation 
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 

	
	// Setting Secretion and Uptake Rates for each metabolites
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

 	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	int glucose_substate_index = microenvironment.find_density_index("glucose");
	int ECM_substate_index = microenvironment.find_density_index("ECM");
	
	
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 0; 

	cell_defaults.phenotype.secretion.uptake_rates[glucose_substate_index] = 0; 
	cell_defaults.phenotype.secretion.secretion_rates[glucose_substate_index] = 1.0; 
	cell_defaults.phenotype.secretion.saturation_densities[glucose_substate_index] = 10.0; 
	
	cell_defaults.phenotype.secretion.uptake_rates[ECM_substate_index] = 0; 
	cell_defaults.phenotype.secretion.secretion_rates[ECM_substate_index] = 1.0; 
	cell_defaults.phenotype.secretion.saturation_densities[ECM_substate_index] = 10.0;   
	
	// Defining Custom Data
	cell_defaults.custom_data.add_variable( "energy", "dimensionless" , 7.01 ); 
	cell_defaults.custom_data.add_variable( "energy_creation_rate", "1/min/mmHg" , 0.01 ); 
	cell_defaults.custom_data.add_variable( "energy_use_rate", "1/min" , 0.02 ); 
	cell_defaults.custom_data.add_variable( "cycle_energy_threshold", "dimensionless" , 10.0 ); 
	cell_defaults.custom_data.add_variable( "death_energy_threshold", "dimensionless" , 7.0 );
	cell_defaults.custom_data.add_variable( "alpha", "none" , 0.0 ); 
	cell_defaults.custom_data.add_variable( "beta", "none" , 0.0 ); 
	cell_defaults.custom_data.add_variable( "gamma", "none" , 0.5 );
	cell_defaults.custom_data.add_variable( "rho", "none" , 0.0 );		
	
 
	// Setting cell volumes, target volumes, and geometry (radii)
/* 	cell_defaults.phenotype.volume.total = 0.0;
	cell_defaults.phenotype.volume.solid = 0.0;
	cell_defaults.phenotype.volume.fluid = 0.0;
	cell_defaults.phenotype.volume.fluid_fraction =	0.0;
	cell_defaults.phenotype.volume.nuclear = 0.0;
	cell_defaults.phenotype.volume.nuclear_fluid = 0.0;
	cell_defaults.phenotype.volume.nuclear_solid = 0.0;
	cell_defaults.phenotype.volume.cytoplasmic = 0.0;
	cell_defaults.phenotype.volume.cytoplasmic_fluid = 0.0;
	cell_defaults.phenotype.volume.cytoplasmic_to_nuclear_ratio = 0.0;
	cell_defaults.phenotype.volume.target_solid_cytoplasmic = 0.0;
	cell_defaults.phenotype.volume.target_solid_nuclear = 0.0;
	cell_defaults.phenotype.volume.target_fluid_fraction = 0.0;
	cell_defaults.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = 0.0;
	cell_defaults.phenotype.volume.cytoplasmic_biomass_change_rate = 0.0;
	cell_defaults.phenotype.volume.nuclear_biomass_change_rate = 0.0;
	cell_defaults.phenotype.volume.fluid_change_rate = 0.0;
	cell_defaults.phenotype.geometry.radius = 0.0;
	cell_defaults.phenotype.geometry.nuclear_radius = 0.0;
	cell_defaults.phenotype.geometry.surface_area = 0.0; */
	 
	
	// Setting Motility
	cell_defaults.phenotype.motility.is_motile = false;
	//cell_defaults.phenotype.motility.migration_speed = 0.0;
	//cell_defaults.phenotype.motility.migration_bias = 0.0;


	// Defining Wound Cell
	wound_cell = cell_defaults; 
	wound_cell.type = 1; 
	wound_cell.name = "wound cell"; 
	wound_cell.parameters.pReference_live_phenotype = &( wound_cell.phenotype );
	wound_cell.phenotype.cycle.data.transition_rate(ncycle_Start_i,ncycle_Start_i) = 0.0;
	wound_cell.functions.update_phenotype = NULL;

	
	
	
	// Setting cell volumes, target volumes, and geometry (radii)
/* 	wound_cell.phenotype.volume.total = 4188.8;
	wound_cell.phenotype.volume.fluid_fraction = 0.7;
	wound_cell.phenotype.volume.fluid = wound_cell.phenotype.volume.fluid_fraction*wound_cell.phenotype.volume.total;
	wound_cell.phenotype.volume.solid =wound_cell.phenotype.volume.total-wound_cell.phenotype.volume.fluid;
	wound_cell.phenotype.volume.cytoplasmic_to_nuclear_ratio = 0.75;
	wound_cell.phenotype.volume.cytoplasmic = 0.75*wound_cell.phenotype.volume.total;
	wound_cell.phenotype.volume.nuclear = wound_cell.phenotype.volume.total-wound_cell.phenotype.volume.cytoplasmic;
	wound_cell.phenotype.volume.nuclear_fluid = 733.04;
	wound_cell.phenotype.volume.nuclear_solid = 314.16;
	wound_cell.phenotype.volume.cytoplasmic = 3141.6;
	wound_cell.phenotype.volume.cytoplasmic_fluid = 2199.12;
	
	wound_cell.phenotype.volume.target_solid_cytoplasmic = 942.48;
	wound_cell.phenotype.volume.target_solid_nuclear = 314.16;
	wound_cell.phenotype.volume.target_fluid_fraction = 0.7;
	wound_cell.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = 0.75;
	wound_cell.phenotype.volume.cytoplasmic_biomass_change_rate = 0.0;
	wound_cell.phenotype.volume.nuclear_biomass_change_rate = 0.0;
	wound_cell.phenotype.volume.fluid_change_rate = 0.0;
	wound_cell.phenotype.geometry.radius = 10;
	wound_cell.phenotype.geometry.nuclear_radius = 6.3;
	wound_cell.phenotype.geometry.surface_area = 1256.637; */
	
	
	// Setting Secretion Rates for Glucose
	wound_cell.phenotype.secretion.secretion_rates[glucose_substate_index]= 0.1; // This should be tuned
	wound_cell.phenotype.secretion.secretion_rates[ECM_substate_index]= 0.0;
	
	
	

	// Defining Bacteria
	bacterial_colony = cell_defaults;
	bacterial_colony.type = 2;
	bacterial_colony.name = "Anaerobic";
	bacterial_colony.parameters.pReference_live_phenotype = &( bacterial_colony2.phenotype );
	bacterial_colony.functions.update_phenotype = energy_update_function;
	bacterial_colony.custom_data["alpha"]=0;
	bacterial_colony.custom_data["beta"]=1;
	bacterial_colony.custom_data["gamma"]=1.0; 
	bacterial_colony.custom_data["rho"]=1.0; //// CHANGE HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// Setting cell volumes, target volumes, and geometry (radii)
/* 	bacterial_colony.phenotype.volume.total = parameters.doubles("bacterial_colony_total_volume");
	bacterial_colony.phenotype.volume.solid = parameters.doubles("bacterial_colony_solid_volume");
	bacterial_colony.phenotype.volume.fluid = parameters.doubles("bacterial_colony_fluid_volume");
	bacterial_colony.phenotype.volume.fluid_fraction = parameters.doubles("bacterial_colony_fluid_fraction");
	bacterial_colony.phenotype.volume.nuclear = parameters.doubles("bacterial_colony_nuclear");
	bacterial_colony.phenotype.volume.nuclear_fluid = parameters.doubles("bacterial_colony_nuclear_fluid");
	bacterial_colony.phenotype.volume.nuclear_solid = parameters.doubles("bacterial_colony_nuclear_solid");
	bacterial_colony.phenotype.volume.cytoplasmic = parameters.doubles("bacterial_colony_cytoplasmic");
	bacterial_colony.phenotype.volume.cytoplasmic_fluid = parameters.doubles("bacterial_colony_cytoplasmic_fluid");
	bacterial_colony.phenotype.volume.cytoplasmic_solid = parameters.doubles("bacterial_colony_cytoplasmic_solid");
	bacterial_colony.phenotype.volume.cytoplasmic_to_nuclear_ratio = parameters.doubles("bacterial_colony_cytoplasmic_to_nuclear_ratio");
	bacterial_colony.phenotype.volume.target_solid_cytoplasmic = parameters.doubles("bacterial_colony_target_solid_cytoplasmic");
	bacterial_colony.phenotype.volume.target_solid_nuclear = parameters.doubles("bacterial_colony_target_solid_nuclear");
	bacterial_colony.phenotype.volume.target_fluid_fraction = parameters.doubles("bacterial_colony_target_fluid_fraction");
	bacterial_colony.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = parameters.doubles("bacterial_colony_target_cytoplasmic_to_nuclear_ratio");
	bacterial_colony.phenotype.volume.cytoplasmic_biomass_change_rate = parameters.doubles("bacterial_colony_cytoplasmic_biomass_change_rate");
	bacterial_colony.phenotype.volume.nuclear_biomass_change_rate = parameters.doubles("bacterial_colony_nuclear_biomass_change_rate");
	bacterial_colony.phenotype.volume.fluid_change_rate = parameters.doubles("bacterial_colony_fluid_change_rate");
	bacterial_colony.phenotype.geometry.radius = parameters.doubles("bacterial_colony_radius");
	bacterial_colony.phenotype.geometry.nuclear_radius = parameters.doubles("bacterial_colony_nuclear_radius");
	bacterial_colony.phenotype.geometry.surface_area = parameters.doubles("bacterial_colony_surface_area");	 */
	bacterial_colony.phenotype.secretion.secretion_rates[ECM_substate_index] = 1.0; 
	bacterial_colony.phenotype.secretion.uptake_rates[glucose_substate_index] =1.0; 	
	//bacterial_colony.functions.volume_update_function = anuclear_volume_model;
	
	
	
	// Defining Bacteria
	bacterial_colony2 = cell_defaults;
	bacterial_colony2.type = 3;
	bacterial_colony2.name = "Aerobic";
	bacterial_colony2.parameters.pReference_live_phenotype = &( bacterial_colony2.phenotype );
	bacterial_colony2.functions.update_phenotype = energy_update_function;

	bacterial_colony2.custom_data["alpha"]=1;
	bacterial_colony2.custom_data["beta"]=0;
	bacterial_colony2.custom_data["gamma"]=1.0; 
	bacterial_colony2.custom_data["rho"]=0.0; 
	// Setting cell volumes, target volumes, and geometry (radii)
/* 	bacterial_colony2.phenotype.volume.total = parameters.doubles("bacterial_colony2_total_volume");
	bacterial_colony2.phenotype.volume.solid = parameters.doubles("bacterial_colony2_solid_volume");
	bacterial_colony2.phenotype.volume.fluid = parameters.doubles("bacterial_colony2_fluid_volume");
	bacterial_colony2.phenotype.volume.fluid_fraction = parameters.doubles("bacterial_colony2_fluid_fraction");
	bacterial_colony2.phenotype.volume.nuclear = parameters.doubles("bacterial_colony2_nuclear");
	bacterial_colony2.phenotype.volume.nuclear_fluid = parameters.doubles("bacterial_colony2_nuclear_fluid");
	bacterial_colony2.phenotype.volume.nuclear_solid = parameters.doubles("bacterial_colony2_nuclear_solid");
	bacterial_colony2.phenotype.volume.cytoplasmic = parameters.doubles("bacterial_colony2_cytoplasmic");
	bacterial_colony2.phenotype.volume.cytoplasmic_fluid = parameters.doubles("bacterial_colony2_cytoplasmic_fluid");
	bacterial_colony2.phenotype.volume.cytoplasmic_solid = parameters.doubles("bacterial_colony2_cytoplasmic_solid");
	bacterial_colony2.phenotype.volume.cytoplasmic_to_nuclear_ratio = parameters.doubles("bacterial_colony2_cytoplasmic_to_nuclear_ratio");
	bacterial_colony2.phenotype.volume.target_solid_cytoplasmic = parameters.doubles("bacterial_colony2_target_solid_cytoplasmic");
	bacterial_colony2.phenotype.volume.target_solid_nuclear = parameters.doubles("bacterial_colony2_target_solid_nuclear");
	bacterial_colony2.phenotype.volume.target_fluid_fraction = parameters.doubles("bacterial_colony2_target_fluid_fraction");
	bacterial_colony2.phenotype.volume.target_cytoplasmic_to_nuclear_ratio = parameters.doubles("bacterial_colony2_target_cytoplasmic_to_nuclear_ratio");
	bacterial_colony2.phenotype.volume.cytoplasmic_biomass_change_rate = parameters.doubles("bacterial_colony2_cytoplasmic_biomass_change_rate");
	bacterial_colony2.phenotype.volume.nuclear_biomass_change_rate = parameters.doubles("bacterial_colony2_nuclear_biomass_change_rate");
	bacterial_colony2.phenotype.volume.fluid_change_rate = parameters.doubles("bacterial_colony2_fluid_change_rate");
	bacterial_colony2.phenotype.geometry.radius = parameters.doubles("bacterial_colony2_radius");
	bacterial_colony2.phenotype.geometry.nuclear_radius = parameters.doubles("bacterial_colony2_nuclear_radius");
	bacterial_colony2.phenotype.geometry.surface_area = parameters.doubles("bacterial_colony2_surface_area");	 */
	bacterial_colony2.phenotype.secretion.secretion_rates[ECM_substate_index] = 0.0;  
	bacterial_colony2.phenotype.secretion.uptake_rates[glucose_substate_index] = 0.5; 	
	bacterial_colony2.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10; 	
	//bacterial_colony2.functions.volume_update_function = anuclear_volume_model;
	
	
	return; 
}




void setup_microenvironment( void )
{
	// set domain parameters 
	
/* now this is in XML 
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
*/
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	// no gradients need for this example 
	microenvironment.add_density( "ECM", "dimensionless", 0.0 , 0.0 ); 	
	microenvironment.add_density( "glucose", "dimensionless", 1.6e3, 0.0 ); 

	default_microenvironment_options.calculate_gradients = true; 
	
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = false;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector_air( 3 ); // 5% o2
	bc_vector_air[0]=20.0;
	bc_vector_air[1]=0.0;
	bc_vector_air[2]=0.0;
	 
	/* std::vector<double> bc_vector_wound( 3 ); // 5% o2
	bc_vector_wound[0]=0.0;
	bc_vector_wound[1]=0.0;
	bc_vector_wound[2]=100; 
	 */
	//default_microenvironment_options.Dirichlet_activation_vector[1] = false;
	//default_microenvironment_options.Dirichlet_activation_vector[1]
	//default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// initialize BioFVM 

initialize_microenvironment(); 	

	
 	for( int n = 0; n < microenvironment.mesh.voxels.size() ; n++ )
	{
		
		microenvironment(n)[1] = 0.0;
		microenvironment(n)[2] = 0.0;
		std::vector<double> position = microenvironment.mesh.voxels[n].center; 
		 if(   position[1] >- 220  )
		{	
		microenvironment.add_dirichlet_node( n,bc_vector_air  );
							
		}
	else
		{
			
		// microenvironment.add_dirichlet_node( n,bc_vector_wound );	
		// microenvironment(n)[2] = 1.0; 
			
		}
		//microenvironment(n)[nECM] = 1.0;  
		
		
		
	}
	
		microenvironment.set_substrate_dirichlet_activation(1,false);
		microenvironment.set_substrate_dirichlet_activation(2,false);
			
	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pC;
	
	// Wound Cell Seeding

	for (int i=-240; i<250; i+=10)
	{
	pC = create_cell(wound_cell); 
	pC->assign_position( i  , -230, 0.0 );
	pC->is_movable=false;
	}

	for (int i=-100; i<100; i+=5)
	{
		
	pC = create_cell(bacterial_colony2); 
	pC->assign_position( i  , -210, 0.0 );
	}
	
	for (int i=-100; i<100; i+=5)
	{
		
	pC = create_cell(bacterial_colony); 
	pC->assign_position( i  , -215, 0.0 );
	}	
/* 	for (int i=-100; i<100; i+=4)
	{
		
	pC = create_cell(bacterial_colony2); 
	pC->assign_position( i  , -200, 0.0 );
	}
	for (int i=-100; i<100; i+=4)
	{
		
	pC = create_cell(bacterial_colony2); 
	pC->assign_position( i  , -205, 0.0 );
	}	
	 */
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
		
	if( pCell->type == 1 )
	{
		 output[0] = "black"; 
		 output[2] = "black"; 
	}
	if( pCell->type == 2)
	{
		 output[0] = "red";
		 output[2] = "red"; 
	}
	
	if( pCell->type == 3)
	{
		 output[0] = "blue";
		 output[2] = "blue"; 
	}	
	
	return output; 
}

/* void anuclear_volume_model( Cell* pCell, Phenotype& phenotype , double dt )
{

	bacterial_colony.phenotype.volume.fluid = bacterial_colony.phenotype.volume.fluid_change_rate * ( bacterial_colony.phenotype.volume.target_fluid_fraction * bacterial_colony.phenotype.volume.total - bacterial_colony.phenotype.volume.fluid) * dt;
	
	bacterial_colony.phenotype.volume.cytoplasmic_solid = bacterial_colony.phenotype.volume.cytoplasmic_biomass_change_rate * (bacterial_colony.phenotype.volume.target_solid_cytoplasmic - bacterial_colony.phenotype.volume.cytoplasmic_solid) * dt;
	
	bacterial_colony.phenotype.volume.cytoplasmic_fluid = bacterial_colony.phenotype.volume.cytoplasmic_fluid;
	bacterial_colony.phenotype.volume.solid = bacterial_colony.phenotype.volume.cytoplasmic_solid;
	bacterial_colony.phenotype.volume.total = bacterial_colony.phenotype.volume.solid + bacterial_colony.phenotype.volume.fluid;
	
	return; 
}
 */


void energy_update_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	static int nO2 = microenvironment.find_density_index( "oxygen" );
	
	static int nGlucose = microenvironment.find_density_index( "glucose" ); 	

	static int nE = pCell->custom_data.find_variable_index( "energy" ); 
	static int nA = pCell->custom_data.find_variable_index( "energy_creation_rate" ); 
	static int nB = pCell->custom_data.find_variable_index( "energy_use_rate" ); 
	
	static int nBirth = pCell->custom_data.find_variable_index( "cycle_energy_threshold" );  
	static int nDeath = pCell->custom_data.find_variable_index( "death_energy_threshold" ); 

	static int nAlpha = pCell->custom_data.find_variable_index( "alpha" ); 
	static int nBeta = pCell->custom_data.find_variable_index( "beta" ); 
	static int nGamma = pCell->custom_data.find_variable_index( "gamma" ); 
	static int nRho = pCell->custom_data.find_variable_index( "rho" );
	
	
	
	double O2 = pCell->nearest_density_vector()[nO2]; 	
	double Glucose = pCell->nearest_density_vector()[nGlucose]; 	
	
	pCell->custom_data[nE] += dt*( pCell->custom_data[nA] * (O2/10) * Glucose * pCell->custom_data[nAlpha] +  pCell->custom_data[nBeta] * pCell->custom_data[nA] * Glucose - pCell->custom_data[nGamma] * pCell->custom_data[nB] - pCell->custom_data[nRho] * (O2/10)); 
	
	//phenotype.secretion.secretion_rates[1] = 1.0; 
//	phenotype.secretion.uptake_rates[2] =0.5; 		
	
	
	static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	static int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	static int ncycle_End_i = live.find_phase_index( PhysiCell_constants::live );
	static int ncycle_Start_i = live.find_phase_index( PhysiCell_constants::live );
	

	// No Proliferation, Apoptosis, and Necrosis
/* 	phenotype.death.rates[necrosis_model_index] = 0.0; 
	phenotype.cycle.data.transition_rate(ncycle_Start_i,ncycle_End_i) = 0.0;
	phenotype.death.rates[apoptosis_model_index] = 0.0; */

	// die if energy is low 
	if( pCell->custom_data[nE] < pCell->custom_data[nDeath] )
	{
		phenotype.death.rates[apoptosis_model_index] = 0.001; 
		
		//return; 
	}

	if( pCell->custom_data[nE] > pCell->custom_data[nBirth])
	{
		phenotype.cycle.data.transition_rate( ncycle_Start_i,ncycle_End_i ) = 0.01; 
		
		//return; 
	}
		
		
		
	return;


	

}





void update_Dirichlet_Nodes(void) 

{
	
	for( int n = 0; n < microenvironment.mesh.voxels.size() ; n++ )
	{
		
		
		if(  microenvironment.nearest_density_vector( n ) [1] > 0 )
		{	
		
		
		microenvironment.remove_dirichlet_node(n);
		
			
		}

		
	}
	
	
} 



 void make_adjustments(void)


{
	
	double dnodes=0;
	double leaked_glucose=0.0;
	for( int n = 0; n < microenvironment.mesh.voxels.size() ; n++ )
	{
		if(microenvironment.is_dirichlet_node(n))
		{ 
		dnodes++;
		
			double  glucose_density =microenvironment.nearest_density_vector(n)[2];//[0];
			//std::cout<< glucose_density;
			if(microenvironment.nearest_density_vector(n)[2]>0.0)
			{
			
			std::vector<double> position = microenvironment.mesh.voxels[n].center;
			//std::cout<< microenvironment.mesh.voxels[n].center;
			double offset=20;
			
			std::vector<std::vector<double>> neighbor_voxels(4);
			neighbor_voxels[0]={position[0]+offset,position[1],position[2]};
			neighbor_voxels[1]={position[0],position[1]+offset,position[2]};
			neighbor_voxels[2]={position[0]-offset,position[1],position[2]};
			neighbor_voxels[3]={position[0],position[1]-offset,position[2]};
			
			double non_air=0.0;
			bool check_it[4]={0,0,0,0};
			for (int m=0;m<4;++m)
				
				{ 
					
					if(fabs(neighbor_voxels[m][0])>240||fabs( neighbor_voxels[m][1])>240)
					{	check_it[m]=0;}
					else 
						
					{ 
						if (!(microenvironment.is_dirichlet_node(microenvironment.nearest_voxel_index(neighbor_voxels[m]))))
					{
						non_air+=1.0;
						check_it[m]=1;

					}
					
					
					}
					
					
				}
			
			for(int j=0;j<4;j++)
	
			{
				if(check_it[j]==1)
				{
				microenvironment(microenvironment.nearest_voxel_index(neighbor_voxels[j]))[2] += (glucose_density/non_air);
				microenvironment(n)[2]-=(glucose_density/non_air);
					//std::cout<<glucose_density/non_air;
					//std::cout<<non_air;
				}
			//std::cout<<non_air;		
			}
			
			if (non_air==0)
			{
					
				leaked_glucose+=microenvironment.nearest_density_vector(n)[2];
				microenvironment.nearest_density_vector(n)[2]=0;
				
			}
		
		  
		//microenvironment(n)[2]=0;
		}// if glucose found
		
			
			
			
		}// end if is_dirichlet_node
		
		
		
		
}// end of for

	double total=microenvironment.mesh.voxels.size();
//std::cout<<leaked_glucose(total-dnodes);
	for( int i = 0; i < microenvironment.mesh.voxels.size() ; i++ )
	{
		
		if (!(microenvironment.is_dirichlet_node(i)))
		{
			
			microenvironment.nearest_density_vector(i)[2]+=(leaked_glucose/(total-dnodes));
			
		}
		
		
		
		
	}
	 

	



/* if( PhysiCell_globals.current_time>60.0)
	{
		
		for( int n = 0; n < microenvironment.mesh.voxels.size() ; n++ )
		{			
			std::vector<double> position = microenvironment.mesh.voxels[n].center;

			if((fabs(position[0])<50&& position[1]<-100&&position[1]>-220))
				{
					
					microenvironment(n)[1] = 1.0;
					
				}
		
		}
	} */

	
}
/* void start_secretion (void)

{
	
	for( int n = 0; n < microenvironment.mesh.voxels.size() ; n++ )
		{			std::vector<double> position = microenvironment.mesh.voxels[n].center;

				if(position[1]<-220)
				{
					
					microenvironment(n)[2] = 1.0;
					
				}
		
		}
	
	
} */