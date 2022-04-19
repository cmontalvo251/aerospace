%%OUTPUT EVERYTHING TO A .MASS and .AERO FILE
MASSFILE = {[num2str(Weight),'	!Weight_(N)_(assume_m1=m2)'];['0	!SLCG(m)'];['0	!BLCG(m)'];['0	!WLCG(m)'];[num2str(Ixx),'        !Ixx(kg*m^2)'];[num2str(Iyy),'	!Iyy'];[num2str(Izz),'	!Izz'];[num2str(II(1,2)),'	!Ixy'];[num2str(II(1,3)),'	!Ixz'];[num2str(II(2,3)),'	!Iyz']}
fid = fopen('UAV.MASS','wb');
for ii = 1:length(MASSFILE)
  fprintf(fid,'%s \n',MASSFILE{ii,:});
end
AEROFILE = {[num2str(CL0_total),'		!C_L_0'];
	    [num2str(CD0_total),' 		!C_D_0'];
	    [num2str(Cm0_total),' 		!C_m_0'];
	    ['0.00000 		!C_D_u'];
	    [num2str(CLA_total),'	        !C_L_alpha'];
	    [num2str(CDA2_total),'     	!C_D_alpha2'];
	    [num2str(Cm_alpha),' 	!C_m_alpha'];
	    ['0.00000 		!C_m_alpha_dot'];
	    ['0.00000 		!C_m_u'];
	    [num2str(CL_q),' 		!C_L_q'];
	    [num2str(Cm_q),'	 	!C_m_q'];
	    [num2str(Clift_delev),' 		!C_L_de'];
	    [num2str(Cm_delev),'     	!C_m_de'];
	    [num2str(C_T),' 		!C_x_delThrust'];
	    [num2str(Cy_beta),' 	!C_y_beta'];
	    [num2str(Cl_beta),'      	!C_l_beta'];
	    [num2str(Cn_beta_total),'      	!C_n_beta'];
	    [num2str(Cl_p),'         !C_l_p'];
	    [num2str(Cn_p),'      	!C_n_p'];
	    [num2str(Cl_r),'		!C_l_r'];
	    [num2str(Cn_r),'		!C_n_r'];
	    ['-0.2559      	!C_l_da'];
	    ['-0.0216		!C_n_da'];
	    [num2str(Cy_drud),'		!C_y_dr'];
	    [num2str(Cl_drud),'      	!C_l_dr'];
	    [num2str(Cn_drud),'		!C_n_dr'];
	    ['0.0057		!C_y_p'];
	    [num2str(S),'		!Reference_Area(m^2)'];
	    [num2str(b),'		!Wingspan(m)'];
	    [num2str(cbar),'		!Mean_chord(m)'];
	    [num2str(Vcruise),'	!Trim_Velocity(m/s)']}

fid = fopen('UAV.AERO','wb');
for ii = 1:length(AEROFILE)
  fprintf(fid,'%s \n',AEROFILE{ii,:});
end