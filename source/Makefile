# Makefile
# erzeugt das Programm FiPPS
# Finite-elemente ProgrammPaket zur Strukturberechnung

include ../makeconfig.mk

PROGRAM := FiPPS
OBJEKTE := \
write_quad_matrix.o \
lsolid20_strain_displ_matrix.o \
lsolid20_eps_th.o \
transform_eg_strain.o \
\
hexa20_shapefunctions_xi_eta_zeta.o \
cosinematrix_esys.o \
quad8s_temploadvec.o \
lsolid20_temploadvec.o \
beam2_shape_functions.o \
beam2_couplingmatrix.o \
\
lsolid20_transformation_matrix.o \
nodalCoordSysTransMat.o \
quad8_normals.o \
quad8_ip_coordsys.o \
quad8_bmat.o \
beam2_cosinematrix.o \
beam2_rotation.o \
hexa20_shapefunctions_x_y_z.o \
transform_eg_cmat.o \
get_node_temperatures.o \
get_beam2_temperatures.o \
get_quad8_temperatures.o \
get_lsolid20_temperatures.o \
vtkoutpointscalar.o \
beam2_temploadvec_control.o \
quad8_temploadvec_control.o \
lsolid20_temploadvec_control.o \
lsolid20_initstressvec.o \
\
mat1_calc_missing_value.o \
mat8.o \
quad8_area.o \
mat20.o \
impc_contact_quad8_pointprojection.o \
impc_contact_lsolid20_naturalcoordinates.o \
quad8_ansatzfunction_xieta.o \
quad8_mpc_factors.o \
spc1_mod_dof.o \
estimate_number_nonzeros_simple.o \
estimate_number_nonzeros_block.o \
assemble_aa.o \
quad8_stiff.o \
beam2_get_properties.o \
beam2_stiff.o \
lsolid20_stiff.o \
blow_up_stiff.o \
load_force.o \
load_moment.o \
load_p2load.o \
load_p8load.o \
load_p20load.o \
temperature_loads.o \
initstress_loads.o \
quad8_geostiff.o \
lsolid20_geostiff.o \
beam2_geostiff.o \
quad8_ip_strain.o \
principalValues.o \
quad8_result_extrapolation_transform.o \
transform_eg_stress.o \
lsolid20_extrapolate_gp_strain.o \
lsolid20_extrapolate_gp_stress.o \
\
boundary_spc_control.o \
boundary_coupling_control.o \
boundary_mpc_control.o \
input_process_line.o \
quad8_prepareproperties.o \
beam2_prepareProperties.o \
lsolid20_prepareProperties.o \
failure_prepare.o \
aeroload_prepare.o \
impc_contact_node_beam2.o \
impc_contact_node_quad8.o \
impc_contact_node_lsolid20.o \
vtkoutcelldata_int.o \
clear_node_node_con.o \
get_node_node_con.o \
estimate_number_nonzeros.o \
quad8_stiff_control.o \
beam2_stiff_control.o \
lsolid20_stiff_control.o \
load_loadvector.o \
vtkoutstatic.o \
quad8_geostiff_control.o \
lsolid20_geostiff_control.o \
beam2_geostiff_control.o \
vtkouteig.o \
vtkoutcelldata.o \
vtkoutcelldata_vector.o \
vtkoutpointvektor.o \
beam2_stress.o \
quad8_results_elem.o \
quad8_results_elem_lay.o \
quad8_results_node.o \
quad8_strainEnergy_control.o \
quad8_strainEnergy.o \
lsolid20_strains_stresses.o \
lsolid20_results_elem_lay.o \
lsolid20_strainEnergy_control.o \
lsolid20_strainEnergy.o \
get_inverse3.o \
get_restored_eigenvektors.o \
\
input_tf.o \
update_geometry.o \
init_values.o \
impc_control.o \
vtkoutprep.o \
vtkoutmesh.o \
vtkoutelementprop.o \
vtkoutelemcoord.o \
vtkoutfixeddispl.o \
boundary_varload_control.o \
init_matrices.o \
init_vectors.o \
init_elem_proc_dist.o \
sol1_get_stiffness_matrices.o \
load_get_loadvector.o \
sol1_solve.o \
sol1_calculate_tse.o \
writeCompTime.o \
sol1_get_utot.o \
sol2_reset_kgeo.o \
sol2_get_geometric_stiffness_matrix.o \
sol2_solve.o \
sol2_extract_eigenvalues.o \
sol2_output_control.o \
sol2_solve_condensed_problem.o \
vtkoutend.o \
free_globals.o \
aero_coupling.o \
load_add_aeropressure.o \
load_add_aeropressure_beam2.o \
load_add_aeropressure_quad8.o \
load_remove_aeropressure.o \
printdisplacements.o \
sol1_output_control.o \
sol1_output_advila.o \
sol1_output_hymowi2.o \
sol1_output_glawi2.o \
sol1_output_glawi.o \
sol1_output_sab.o \
sol1_reaction_forces.o \
\
main.o \

# koopt_sol1_output.o \
# mat_transform_mat_shell.o \
# sol1_output_glawi.o \
# sol1_output_hymowi.o \
# sol1_output_koopt.o \
# sol1_output_optitube.o \
# sol1_output_sab.o \
# writeinput.o \
# write_matrix.o \
# write_quad_matrix.o

###############################################################################
## rules follow
###############################################################################
.PHONY: default clean

all: $(PROGRAM)

$(PROGRAM): kompilieren move

kompilieren: 
	-${FLINKER} -o $(PROGRAM) $(FIPPS2_LIB) $(APAME_LIB) $(PANEL2D_LIB) $(XFOILWRAPPER_LIB) ${SLEPC_LIB} -L../include/ -lmodule -L../ -luseroutput

move:
	$(MV) $(PROGRAM) ../.

#flib: $(OBJEKTE)

flib: $(OBJEKTE)
	$(AR) $(AR_FLAGS) $(FIPPS2_LIB) $^

clean2: clean
	$(RM) $(FIPPS2_LIB)
