function [cikan] = fn_silver_k_2(lamda)
lamda=lamda*10^6;
data=[2.700E-01	1.318E+00
2.800E-01	1.240E+00
2.900E-01	1.168E+00
3.000E-01	9.126E-01
3.100E-01	5.192E-01
3.200E-01	3.334E-01
3.300E-01	6.779E-01
3.400E-01	1.084E+00
3.500E-01	1.341E+00
3.600E-01	1.536E+00
3.699E-01	1.699E+00
3.799E-01	1.841E+00
3.899E-01	1.971E+00
3.999E-01	2.090E+00
4.099E-01	2.203E+00
4.199E-01	2.310E+00
4.299E-01	2.413E+00
4.399E-01	2.513E+00
4.499E-01	2.609E+00
4.599E-01	2.703E+00
4.699E-01	2.795E+00
4.799E-01	2.886E+00
4.899E-01	2.974E+00
4.999E-01	3.062E+00
5.099E-01	3.148E+00
5.199E-01	3.233E+00
5.299E-01	3.318E+00
5.399E-01	3.401E+00
5.499E-01	3.484E+00
5.599E-01	3.566E+00
5.699E-01	3.648E+00
5.799E-01	3.728E+00
5.899E-01	3.809E+00
5.999E-01	3.889E+00
6.099E-01	3.968E+00
6.199E-01	4.048E+00
6.299E-01	4.126E+00
6.399E-01	4.205E+00
6.499E-01	4.283E+00
6.599E-01	4.361E+00
6.699E-01	4.438E+00
6.799E-01	4.516E+00
6.899E-01	4.593E+00
6.999E-01	4.669E+00
7.099E-01	4.746E+00
7.199E-01	4.823E+00
7.299E-01	4.899E+00
7.399E-01	4.975E+00
7.499E-01	5.051E+00
7.599E-01	5.126E+00
7.699E-01	5.202E+00
7.799E-01	5.278E+00
7.899E-01	5.353E+00
7.999E-01	5.428E+00
8.099E-01	5.503E+00
8.199E-01	5.578E+00
8.299E-01	5.653E+00
8.399E-01	5.728E+00
8.499E-01	5.802E+00
8.599E-01	5.877E+00
8.699E-01	5.951E+00
8.799E-01	6.026E+00
8.899E-01	6.100E+00
8.999E-01	6.174E+00
9.099E-01	6.248E+00
9.199E-01	6.322E+00
9.299E-01	6.396E+00
9.399E-01	6.470E+00
9.499E-01	6.544E+00
9.599E-01	6.618E+00
9.698E-01	6.691E+00
9.799E-01	6.765E+00
9.899E-01	6.838E+00
9.999E-01	6.912E+00
1.010E+00	6.985E+00
1.020E+00	7.059E+00
1.030E+00	7.132E+00
1.040E+00	7.206E+00
1.050E+00	7.279E+00
1.060E+00	7.352E+00
1.070E+00	7.425E+00
1.080E+00	7.498E+00
1.090E+00	7.571E+00
1.100E+00	7.645E+00
1.110E+00	7.718E+00
1.120E+00	7.791E+00
1.130E+00	7.863E+00
1.140E+00	7.936E+00
1.150E+00	8.009E+00
1.160E+00	8.082E+00
1.170E+00	8.155E+00
1.180E+00	8.228E+00
1.190E+00	8.300E+00
1.200E+00	8.373E+00
1.210E+00	8.446E+00
1.220E+00	8.518E+00
1.230E+00	8.591E+00
1.240E+00	8.664E+00
1.250E+00	8.736E+00
1.260E+00	8.809E+00
1.268E+00	8.869E+00
1.270E+00	8.881E+00
1.278E+00	8.942E+00
1.280E+00	8.954E+00
1.288E+00	9.015E+00
1.290E+00	9.026E+00
1.299E+00	9.090E+00
1.300E+00	9.099E+00
1.309E+00	9.166E+00
1.310E+00	9.171E+00
1.320E+00	9.243E+00
1.328E+00	9.302E+00
1.330E+00	9.316E+00
1.339E+00	9.381E+00
1.340E+00	9.388E+00
1.347E+00	9.442E+00
1.350E+00	9.460E+00
1.359E+00	9.524E+00
1.367E+00	9.586E+00
1.379E+00	9.670E+00
1.388E+00	9.734E+00
1.400E+00	9.820E+00
1.409E+00	9.886E+00
1.418E+00	9.953E+00
1.427E+00	1.002E+01
1.430E+00	1.004E+01
1.437E+00	1.009E+01
1.440E+00	1.011E+01
1.450E+00	1.018E+01
1.460E+00	1.025E+01
1.469E+00	1.032E+01
1.470E+00	1.033E+01
1.480E+00	1.040E+01
1.490E+00	1.047E+01
1.497E+00	1.052E+01
1.500E+00	1.054E+01
1.507E+00	1.060E+01
1.510E+00	1.062E+01
1.518E+00	1.067E+01
1.520E+00	1.069E+01
1.528E+00	1.075E+01
1.530E+00	1.076E+01
1.539E+00	1.083E+01
1.540E+00	1.083E+01
1.547E+00	1.088E+01
1.550E+00	1.090E+01
1.558E+00	1.096E+01
1.560E+00	1.098E+01
1.569E+00	1.104E+01
1.570E+00	1.105E+01
1.577E+00	1.110E+01
1.580E+00	1.112E+01
1.588E+00	1.118E+01
1.590E+00	1.119E+01
1.596E+00	1.124E+01
1.600E+00	1.126E+01
1.608E+00	1.132E+01
1.610E+00	1.134E+01
1.616E+00	1.138E+01
1.620E+00	1.141E+01
1.628E+00	1.147E+01
1.630E+00	1.148E+01
1.636E+00	1.153E+01
1.640E+00	1.155E+01
1.649E+00	1.162E+01
1.650E+00	1.162E+01
1.657E+00	1.168E+01
1.660E+00	1.169E+01
1.666E+00	1.174E+01
1.670E+00	1.177E+01
1.679E+00	1.183E+01
1.680E+00	1.184E+01
1.688E+00	1.189E+01
1.690E+00	1.191E+01
1.696E+00	1.196E+01
1.700E+00	1.198E+01
1.705E+00	1.202E+01
1.710E+00	1.205E+01
1.719E+00	1.212E+01
1.720E+00	1.213E+01
1.728E+00	1.218E+01
1.730E+00	1.220E+01
1.737E+00	1.225E+01
1.740E+00	1.227E+01
1.747E+00	1.232E+01
1.750E+00	1.234E+01
1.756E+00	1.239E+01
1.760E+00	1.241E+01
1.766E+00	1.246E+01
1.770E+00	1.248E+01
1.775E+00	1.252E+01
1.780E+00	1.256E+01
1.785E+00	1.259E+01
1.790E+00	1.263E+01
1.795E+00	1.267E+01
1.800E+00	1.270E+01
1.805E+00	1.274E+01
1.810E+00	1.277E+01
1.815E+00	1.281E+01
1.820E+00	1.284E+01
1.825E+00	1.288E+01
1.830E+00	1.291E+01
1.836E+00	1.296E+01
1.840E+00	1.299E+01
1.846E+00	1.303E+01
1.850E+00	1.306E+01
1.857E+00	1.311E+01
1.860E+00	1.313E+01
1.868E+00	1.318E+01
1.870E+00	1.320E+01
1.878E+00	1.326E+01
1.880E+00	1.327E+01
1.889E+00	1.334E+01
1.890E+00	1.334E+01
1.895E+00	1.338E+01
1.900E+00	1.342E+01
1.906E+00	1.346E+01
1.910E+00	1.349E+01
1.917E+00	1.354E+01
1.920E+00	1.356E+01
1.929E+00	1.362E+01
1.930E+00	1.363E+01
1.934E+00	1.366E+01
1.940E+00	1.370E+01
1.946E+00	1.375E+01
1.950E+00	1.377E+01
1.958E+00	1.383E+01
1.960E+00	1.384E+01
1.970E+00	1.392E+01
1.976E+00	1.396E+01
1.980E+00	1.399E+01
1.988E+00	1.405E+01
1.990E+00	1.406E+01
1.994E+00	1.409E+01
2.000E+00	1.413E+01
2.006E+00	1.418E+01
2.019E+00	1.427E+01
2.025E+00	1.431E+01
2.038E+00	1.440E+01
2.044E+00	1.445E+01
2.057E+00	1.454E+01
2.064E+00	1.459E+01
2.077E+00	1.468E+01
2.084E+00	1.473E+01
2.097E+00	1.483E+01
2.104E+00	1.488E+01
2.118E+00	1.497E+01
2.125E+00	1.502E+01
2.139E+00	1.512E+01
2.146E+00	1.517E+01
2.153E+00	1.523E+01
2.167E+00	1.533E+01
2.175E+00	1.538E+01
2.189E+00	1.548E+01
2.197E+00	1.554E+01
2.204E+00	1.559E+01
2.219E+00	1.570E+01
2.227E+00	1.575E+01
2.235E+00	1.581E+01
2.242E+00	1.586E+01
2.258E+00	1.597E+01
2.266E+00	1.603E+01
2.274E+00	1.609E+01
2.282E+00	1.614E+01
2.298E+00	1.626E+01
2.306E+00	1.632E+01
2.314E+00	1.638E+01
2.323E+00	1.644E+01
2.339E+00	1.656E+01
2.348E+00	1.662E+01
2.356E+00	1.668E+01
2.365E+00	1.674E+01
2.374E+00	1.680E+01
2.382E+00	1.686E+01
2.391E+00	1.692E+01
2.409E+00	1.705E+01
2.418E+00	1.712E+01
2.427E+00	1.718E+01
2.436E+00	1.724E+01
2.445E+00	1.731E+01
2.455E+00	1.738E+01
2.464E+00	1.744E+01
2.473E+00	1.751E+01
2.483E+00	1.758E+01
2.492E+00	1.765E+01
2.502E+00	1.771E+01
2.512E+00	1.778E+01
2.522E+00	1.785E+01
2.531E+00	1.792E+01
2.541E+00	1.799E+01
2.551E+00	1.806E+01
2.561E+00	1.814E+01
2.572E+00	1.821E+01
2.582E+00	1.828E+01
2.592E+00	1.835E+01
2.603E+00	1.843E+01
2.613E+00	1.850E+01
2.624E+00	1.858E+01
2.634E+00	1.865E+01
2.645E+00	1.873E+01
2.656E+00	1.881E+01
2.667E+00	1.889E+01
2.678E+00	1.896E+01
2.689E+00	1.904E+01
2.700E+00	1.912E+01
2.711E+00	1.920E+01
2.723E+00	1.928E+01
2.734E+00	1.936E+01
2.746E+00	1.945E+01
2.758E+00	1.953E+01
2.769E+00	1.961E+01
2.781E+00	1.970E+01
2.793E+00	1.978E+01
2.805E+00	1.987E+01
2.817E+00	1.996E+01
2.830E+00	2.004E+01
2.842E+00	2.013E+01
2.855E+00	2.022E+01
2.867E+00	2.031E+01
2.880E+00	2.040E+01
2.893E+00	2.049E+01
2.906E+00	2.058E+01
2.919E+00	2.068E+01
2.932E+00	2.077E+01
2.946E+00	2.086E+01
2.959E+00	2.096E+01
2.973E+00	2.106E+01
2.986E+00	2.115E+01
3.000E+00	2.125E+01
3.014E+00	2.135E+01
3.028E+00	2.145E+01
3.042E+00	2.155E+01
3.057E+00	2.165E+01
3.071E+00	2.175E+01
3.086E+00	2.186E+01
3.101E+00	2.196E+01
3.116E+00	2.207E+01
3.131E+00	2.217E+01
3.146E+00	2.228E+01
3.161E+00	2.239E+01
3.177E+00	2.250E+01
3.192E+00	2.261E+01
3.208E+00	2.272E+01
3.224E+00	2.284E+01
3.240E+00	2.295E+01
3.256E+00	2.307E+01
3.273E+00	2.318E+01
3.289E+00	2.330E+01
3.306E+00	2.342E+01
3.323E+00	2.354E+01
3.340E+00	2.366E+01
3.358E+00	2.378E+01
3.375E+00	2.390E+01
3.393E+00	2.403E+01
3.411E+00	2.416E+01
3.429E+00	2.428E+01
3.447E+00	2.441E+01
3.465E+00	2.454E+01
3.484E+00	2.467E+01
3.503E+00	2.481E+01
3.522E+00	2.494E+01
3.541E+00	2.508E+01
3.561E+00	2.521E+01
3.580E+00	2.535E+01
3.600E+00	2.549E+01
3.620E+00	2.563E+01
3.641E+00	2.578E+01
3.661E+00	2.592E+01
3.682E+00	2.607E+01
3.703E+00	2.622E+01
3.724E+00	2.637E+01
3.746E+00	2.652E+01
3.768E+00	2.667E+01
3.790E+00	2.683E+01
3.812E+00	2.698E+01
3.834E+00	2.714E+01
3.857E+00	2.730E+01
3.880E+00	2.747E+01
3.904E+00	2.763E+01
3.927E+00	2.780E+01
3.951E+00	2.797E+01
3.976E+00	2.814E+01
4.000E+00	2.831E+01
4.025E+00	2.848E+01
4.050E+00	2.866E+01
4.076E+00	2.884E+01
4.101E+00	2.902E+01
4.128E+00	2.920E+01
4.154E+00	2.939E+01
4.181E+00	2.958E+01
4.208E+00	2.977E+01
4.235E+00	2.996E+01
4.263E+00	3.016E+01
4.292E+00	3.035E+01
4.320E+00	3.056E+01
4.349E+00	3.076E+01
4.379E+00	3.096E+01
4.408E+00	3.117E+01
4.439E+00	3.138E+01
4.469E+00	3.160E+01
4.500E+00	3.182E+01
4.532E+00	3.204E+01
4.564E+00	3.226E+01
4.596E+00	3.249E+01
4.629E+00	3.272E+01
4.662E+00	3.295E+01
4.696E+00	3.318E+01
4.730E+00	3.342E+01
4.765E+00	3.367E+01
4.800E+00	3.391E+01
4.836E+00	3.416E+01
4.872E+00	3.442E+01
4.909E+00	3.467E+01
4.947E+00	3.493E+01
4.985E+00	3.520E+01
5.023E+00	3.547E+01
5.063E+00	3.574E+01
5.103E+00	3.602E+01
5.143E+00	3.630E+01
5.184E+00	3.659E+01
5.226E+00	3.688E+01
5.269E+00	3.717E+01
5.312E+00	3.747E+01
5.356E+00	3.778E+01
5.400E+00	3.809E+01
5.446E+00	3.840E+01
5.492E+00	3.872E+01
5.539E+00	3.905E+01
5.586E+00	3.938E+01
5.635E+00	3.971E+01
5.684E+00	4.005E+01
5.735E+00	4.040E+01
5.786E+00	4.076E+01
5.838E+00	4.112E+01
5.891E+00	4.148E+01
5.945E+00	4.185E+01
6.000E+00	4.223E+01
6.056E+00	4.262E+01
6.114E+00	4.301E+01
6.172E+00	4.341E+01
6.231E+00	4.382E+01
6.292E+00	4.424E+01
6.353E+00	4.466E+01
6.416E+00	4.509E+01
6.480E+00	4.553E+01
6.546E+00	4.598E+01
6.612E+00	4.644E+01
6.681E+00	4.690E+01
6.750E+00	4.738E+01
6.821E+00	4.786E+01
6.894E+00	4.836E+01
6.968E+00	4.886E+01
7.044E+00	4.938E+01
7.121E+00	4.990E+01
7.200E+00	5.044E+01
7.281E+00	5.099E+01
7.364E+00	5.155E+01
7.449E+00	5.212E+01
7.535E+00	5.271E+01
7.624E+00	5.331E+01
7.715E+00	5.392E+01
7.808E+00	5.454E+01
7.903E+00	5.518E+01
8.001E+00	5.584E+01
8.100E+00	5.651E+01
8.203E+00	5.720E+01
8.308E+00	5.790E+01
8.416E+00	5.862E+01
8.527E+00	5.936E+01
8.641E+00	6.011E+01
8.757E+00	6.089E+01
8.877E+00	6.168E+01
9.001E+00	6.250E+01
9.127E+00	6.333E+01
9.257E+00	6.419E+01
9.392E+00	6.507E+01
9.530E+00	6.598E+01
9.672E+00	6.691E+01
9.819E+00	6.787E+01
9.970E+00	6.885E+01
1.013E+01	6.986E+01
1.029E+01	7.090E+01
1.045E+01	7.197E+01
1.062E+01	7.307E+01
1.080E+01	7.421E+01
1.098E+01	7.538E+01
1.117E+01	7.658E+01
1.137E+01	7.783E+01
1.157E+01	7.911E+01
1.178E+01	8.043E+01
1.200E+01	8.180E+01
1.223E+01	8.321E+01
1.246E+01	8.467E+01
1.271E+01	8.618E+01
1.296E+01	8.774E+01
1.323E+01	8.936E+01
1.350E+01	9.103E+01
1.379E+01	9.277E+01
1.409E+01	9.457E+01
1.440E+01	9.643E+01
1.473E+01	9.837E+01
1.507E+01	1.004E+02
1.543E+01	1.025E+02
1.581E+01	1.046E+02
1.620E+01	1.069E+02
1.662E+01	1.093E+02
1.705E+01	1.117E+02
1.751E+01	1.143E+02
1.800E+01	1.169E+02
1.852E+01	1.197E+02
1.906E+01	1.226E+02
1.964E+01	1.256E+02
2.025E+01	1.288E+02
2.090E+01	1.321E+02
2.160E+01	1.356E+02
2.235E+01	1.393E+02
2.314E+01	1.431E+02
2.400E+01	1.471E+02
2.502E+01	1.514E+02];

cikan=interp1(data(:,1),data(:,2),lamda);
end