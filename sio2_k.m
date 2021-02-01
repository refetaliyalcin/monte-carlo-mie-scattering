function result=sio2_k(x)
x=10^6*x;
data=[3.00653594771242e-002	7.12929515454867e-002
3.79084967320261e-002	1.33489783491451e-001
4.96732026143791e-002	2.85231268014320e-001
5.62091503267974e-002	4.38103841136116e-001
7.45098039215686e-002	7.67899985694281e-001
9.67320261437909e-002	7.67899985694281e-001
1.13725490196078e-001	7.67899985694281e-001
1.21568627450980e-001	5.89670388029277e-001
1.29411764705882e-001	6.24739441108538e-002
1.42483660130719e-001	8.90883239609268e-003
1.50326797385621e-001	9.13215925902242e-004
1.54248366013072e-001	1.14116296898087e-004
1.58169934640523e-001	2.05034978483135e-005
1.63398692810458e-001	4.34503052210942e-006
1.67320261437908e-001	1.51081563140158e-006
1.88235294117647e-001	1.16015530173997e-006
2.07843137254902e-001	1.16015530173997e-006
2.24836601307190e-001	5.42958985962062e-007
2.41830065359477e-001	9.83629491042605e-007
2.69281045751634e-001	6.35136956341499e-008
3.09803921568627e-001	7.01258496436122e-008
3.59477124183007e-001	1.35711448982643e-007
3.80392156862745e-001	1.04212819732322e-007
3.98692810457516e-001	6.78486512204820e-008
4.81045751633987e-001	7.74263682681127e-008
5.62091503267974e-001	7.74263682681127e-008
6.22222222222222e-001	7.49120999259123e-008
6.58823529411765e-001	1.07710505603677e-007
6.79738562091503e-001	9.43866113179867e-008
8.03921568627451e-001	1.22915240365122e-007
8.43137254901961e-001	1.11325584003956e-007
8.81045751633987e-001	4.87721659688546e-008
9.22875816993464e-001	9.43866113179867e-008
9.60784313725490e-001	9.75545010046874e-008
1.00000000000000e+000	8.54869143640455e-008
1.21932114882507e+000	4.00083404924448e-008
1.42036553524804e+000	7.74263682681127e-008
1.56657963446475e+000	4.71883862356841e-008
1.93211488250653e+000	4.16938197552849e-007
2.57180156657963e+000	2.32055249753435e-006
3.32114882506527e+000	4.49086240639031e-006
3.79634464751958e+000	4.83713060388846e-005
4.47258485639687e+000	3.07224040918003e-004
4.76501305483029e+000	4.71883862356840e-004
5.00261096605744e+000	5.08268494006711e-003
6.48302872062663e+000	6.61894131365520e-003
7.10443864229765e+000	1.01664296272781e-002
7.48825065274151e+000	2.47893553464351e-002
7.67101827676240e+000	6.45707517032943e-002
7.92689295039165e+000	1.37970089569798e-001
9.02349869451697e+000	2.69220128298041e+000
9.75456919060052e+000	7.18837664294370e-001
1.01018276762402e+001	3.04698957090351e-001
1.08694516971279e+001	1.85702266481110e-001
1.13080939947781e+001	1.52333597657718e-001
1.20391644908616e+001	2.75968946102348e-001
1.26057441253264e+001	3.71443252152732e-001
1.31540469973890e+001	2.26380340952145e-001
1.37937336814621e+001	7.87149724293484e-002
1.41775456919060e+001	4.79737407889900e-002
1.44699738903394e+001	4.34503052210942e-002
1.50000000000000e+001	5.47458564718823e-002
1.56703470031546e+001	8.05953065382090e-002
1.62066246056782e+001	1.22988967837561e-001
1.67429022082019e+001	1.68862595693521e-001
1.74132492113565e+001	2.19915689349749e-001
1.84858044164038e+001	3.05956480048725e-001
1.94242902208202e+001	4.66891600772887e-001
1.99605678233438e+001	7.91883999516008e-001
2.03627760252366e+001	1.36095326627570e+000
2.08990536277603e+001	1.86857492343883e+000
2.14353312302839e+001	2.33897110450336e+000
2.25078864353312e+001	1.32547315675008e+000
2.31782334384858e+001	8.91839712197086e-001
2.38485804416404e+001	5.84427174658512e-001
2.46529968454259e+001	3.26843870030386e-001
2.51892744479495e+001	2.41217703139677e-001
2.55914826498423e+001	2.08598731667757e-001
2.61277602523659e+001	1.90177309853725e-001
2.74684542586751e+001	1.73382689789336e-001
3.00157728706625e+001	1.44111888304600e-001
3.84621451104101e+001	6.78797074257425e-002
4.98580441640379e+001	2.72864139111488e-002];
result=interp1(data(:,1),data(:,2),x);