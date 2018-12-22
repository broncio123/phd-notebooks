#include <pdb.hpp>	
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <algorithm>
#include <iomanip>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/count.hpp>

using namespace std;
using namespace boost;
using namespace boost::accumulators;

#define _USE_MATH_DEFINES

double mempot[20][34] = {
	{2.495,2.5773,2.6103,2.6452,2.4797,2.4579,2.5666,2.3904,2.2742,2.243,2.1455,2.1863,2.0682,2.0939,2.1043,1.9889,2.1038,2.1392,2.0358,2.0304,2.0944,2.1321,2.1932,2.1324,2.28,2.1965,2.2755,2.3471,2.7097,2.4487,2.537,2.5084,2.6783,2.5451},
	{2.6123,2.3606,2.5193,2.425,2.3095,2.3535,2.7121,2.8041,3.2166,3.5115,3.9505,5.2032,4.8477,5.6912,6.0426,5.9402,5.9956,5.135,5.2874,5.425,5.0775,5.0499,4.3479,4.725,4.079,3.8059,3.5937,3.8635,3.5238,3.0055,3.1864,3.607,3.4475,3.0692},
	{3.2655,3.0935,3.1159,3.2158,3.1032,3.151,3.3467,3.7522,4.0959,4.2619,3.9505,3.8533,3.9094,3.9865,4.1968,4.0431,4.6093,4.3083,4.5142,4.5495,4.5179,4.3567,3.6865,3.9874,3.718,4.0821,3.5424,3.2275,2.8745,3.1761,2.7344,2.8674,3.1338,3.0909},
	{2.8153,3.0155,2.8022,2.9281,3.1424,3.7308,4.2977,3.8903,4.789,4.5301,5.2234,4.7512,5.099,4.6796,5.0618,5.4293,5.7079,5.9823,4.6813,4.9549,5.211,4.8267,4.8179,5.177,4.3241,4.0821,4.0958,3.6628,3.0538,3.0872,3.004,2.7395,2.7543,2.7544},
	{5.0691,6.0113,4.8075,5.7807,5.9558,5.4536,4.9337,4.6212,4.1424,4.7925,4.1736,4.2588,4.56,4.8803,4.0967,3.8607,4.3862,4.2477,4.0347,4.0387,4.1124,4.1849,4.9232,4.4148,5.1351,5.7869,5.0341,5.2134,5.2874,4.8579,5.3554,5.2653,4.3125,4.3757},
	{3.246,2.9667,3.4576,3.5121,3.6205,3.5077,3.6344,3.928,3.8196,4.0505,4.6844,4.9519,5.099,4.7749,5.0618,5.0929,4.6963,4.7783,4.6813,4.7318,5.211,4.7314,3.9016,3.9448,3.8472,3.9951,3.4246,3.4957,3.5528,3.1087,3.1041,3.1252,2.781,3.0506},
	{2.4198,2.5455,2.8022,2.9875,3.2703,3.5077,3.4933,3.8903,4.1912,4.0993,4.6044,4.5842,4.7423,5.6912,4.7433,4.9593,4.6963,5.4714,4.9997,4.2012,5.0775,4.2961,4.176,4.2897,3.9956,4.0821,3.8301,3.6628,3.4416,3.2997,3.0281,2.938,2.4407,2.6611},
	{2.6062,2.7726,2.7104,2.6167,2.6236,2.5413,2.7364,2.8554,2.633,2.5842,2.3653,2.4306,2.3444,2.2255,2.305,2.4237,2.274,2.2288,2.1964,2.2386,2.2666,2.3586,2.4291,2.4133,2.5798,2.4458,2.4691,2.2466,2.1163,2.2188,2.1466,2.3475,2.1724,2.5185},
	{3.7505,3.814,3.6624,3.8836,3.7586,3.6958,3.8351,4.0507,4.0088,4.387,4.3361,4.4411,4.647,4.4384,4.8387,4.9593,5.3025,5.135,4.5943,4.3953,4.3237,4.1336,4.3479,4.0783,4.079,3.9951,4.0958,3.6329,3.6134,3.6792,3.4836,3.7071,3.6193,4.1545},
	{2.9041,2.9667,3.2271,3.1417,3.0655,2.7192,2.6536,2.504,2.3506,2.3501,2.1721,2.2893,2.0682,2.1359,2.0978,2.1712,2.1743,2.0837,2.0888,2.1071,2.12,2.301,2.3164,2.2325,2.3167,2.6514,2.5966,2.682,2.8895,3.2484,3.0781,3.3935,3.3009,2.9771},
	{2.5471,2.6614,2.3108,2.4367,2.3449,2.3357,1.9834,2.1427,1.8786,1.864,1.8861,1.6357,1.7317,1.6884,1.7075,1.8184,1.6475,1.8028,1.6993,1.7874,1.7232,1.831,1.9846,1.8654,1.8579,1.9367,2.0495,2.1621,2.317,2.6198,2.6146,2.6813,2.7284,2.5889},
	{2.6975,2.3869,2.4494,2.4849,2.66,3.0026,2.7489,3.0273,3.6316,3.7992,4.6844,5.3574,5.9463,5.6912,4.944,5.247,5.4848,5.9823,7.0792,5.0885,4.9597,6.4362,5.1745,5.177,4.3813,4.1775,4.0958,3.8635,3.4957,3.3264,3.5636,3.3194,3.827,2.9845},
	{3.8867,3.8712,3.9202,4.1068,3.6533,3.7308,3.5756,3.2943,3.2951,3.2664,3.2573,3.3649,3.3072,3.2489,3.5303,3.1268,3.2024,3.1296,3.5528,3.2732,3.4193,3.3226,3.1692,3.4093,3.3859,3.2612,3.4936,3.5484,3.5827,3.3538,4.0337,3.4735,3.907,3.8034},
	{3.4692,3.7087,3.6179,3.2418,3.1032,2.8289,2.6423,2.4588,2.5522,2.4222,2.4332,2.4217,2.3815,2.5449,2.6086,2.4744,2.4595,2.3535,2.7485,2.5917,2.3948,2.3253,2.3002,2.4133,2.4871,2.3774,2.326,2.682,2.6133,2.93,3.2452,2.9139,3.173,3.2738},
	{3.0107,2.7532,2.8822,3.0082,3.2932,3.1118,3.3467,3.2741,3.3157,3.6611,3.8735,3.7479,3.4896,3.4399,3.6755,3.8199,3.7984,3.5255,3.678,3.8989,3.8248,3.6029,3.6547,3.3162,3.1892,3.2414,2.9546,2.8966,2.9683,2.7942,2.5222,2.6085,2.5853,2.7505},
	{2.7617,2.9202,2.9468,2.3685,2.7106,2.8435,2.9565,2.8957,3.076,3.0697,2.995,3.0887,3.0376,2.8878,2.9669,2.8193,2.8747,3.2307,3.0719,2.8003,3.0138,2.8526,2.8863,3.0285,3.0401,3.0301,2.9138,2.8552,2.7751,2.7477,2.552,2.4927,2.703,2.8816},
	{2.929,2.9202,2.6591,2.9673,2.8956,2.9687,2.9878,3.1431,2.9972,2.8904,3.0422,3.0383,3.0019,2.9186,2.9071,2.8491,2.8316,3.0204,2.936,3.0641,3.0965,2.9396,2.9308,2.9957,2.8707,2.9247,3.0416,2.7773,2.9203,2.8102,2.9805,2.9139,3.0597,2.8212},
	{5.2514,5.095,4.5562,4.1713,3.5887,3.4255,3.2808,3.1257,3.4257,3.4315,3.3407,3.4856,3.9539,3.8994,3.5859,4.0943,4.261,3.862,3.7833,4.2618,3.8611,3.6029,3.4829,3.2727,3.2072,3.0623,3.2059,3.3526,3.3656,3.6051,4.1026,3.9435,4.3125,4.2117},
	{3.858,3.4086,3.1981,3.4782,3.0842,3.1712,3.3023,2.9376,3.076,3.0175,3.614,3.4356,3.4614,4.0818,3.8831,3.8607,4.0497,3.9454,3.747,3.7763,3.3068,3.5184,3.1886,3.4093,2.9804,3.0301,2.9686,3.0818,3.2725,3.0872,3.3077,3.1252,3.0597,3.2771},
	{2.6561,3.0155,3.0646,2.8189,2.994,2.8004,2.5983,2.581,2.4247,2.2993,2.1587,2.1452,2.3536,2.2653,2.1108,2.1043,2.1384,2.031,2.1375,2.1746,2.2443,2.1595,2.2301,2.2946,2.4783,2.3857,2.688,2.7049,2.7887,2.8264,2.935,2.988,2.9262,2.6434}};

double betamempot[20][34] = {
	{2.1883,2.5123,2.5123,2.3191,3.434,2.5813,2.8134,2.6194,2.9191,2.9277,2.6271,2.7651,2.6527,2.4372,2.2463,1.9134,1.8433,2.1743,2.1691,2.1367,2.5762,2.3514,2.7033,2.6198,2.9704,3.054,2.6415,2.8434,2.8674,2.6794,2.5353,2.8094,2.6696,2.5063},
	{3.3277,3.6109,3.6109,4.1109,2.7408,5.2204,3.1236,2.8772,3.0732,3.5695,3.3202,2.7651,3.4259,3.1304,2.7651,3.1406,3.042,3.0704,3.6424,3.399,3.2122,2.8702,2.991,3.313,2.6603,2.5532,2.5123,2.4235,2.7804,2.6053,2.5704,2.3711,2.5863,3.2315},
	{3.4613,3.6109,3.6109,2.5014,3.2108,2.2246,2.1893,2.2447,2.5296,2.8276,3.9488,3.3449,3.1382,2.9678,3.3449,3.3348,3.2139,3.7635,3.4012,3.0625,3.2693,3.2558,3.2023,3.4561,2.9704,2.4944,2.3453,2.29,2.9627,2.6053,2.4369,2.7606,2.4378,2.563},
	{2.5168,2.2246,1.665,1.713,1.6848,1.5568,1.8579,2.4146,2.6314,3.0995,2.9838,3.2759,3.3458,3.3245,3.419,2.6416,3.5886,3.2935,3.5553,3.5421,3.1581,2.911,3.2023,3.382,3.2282,2.7916,2.3199,2.368,2.3475,2.2487,2.5353,2.1416,2.2788,2.4018},
	{5.4072,3.6109,3.6109,4.1109,4.8203,5.2204,5.5215,5.7104,5.9636,5.8721,6.0283,5.9839,5.9108,5.9636,5.9839,5.9738,5.9865,6.0661,6.0403,6.107,6.1026,6.089,6.0355,6.021,5.3683,6.0497,6.0088,5.9789,5.9584,5.9375,5.9026,5.8051,5.8051,5.5588},
	{3.2099,2.9178,3.6109,3.4177,3.2108,2.5123,3.0366,3.4078,3.3245,3.1641,3.6304,3.419,3.7136,3.3245,3.499,4.1821,3.6839,3.5012,3.4753,3.399,3.33,3.6911,3.3964,3.0253,2.7656,2.6157,2.7899,2.8434,3.1858,3.2295,2.9069,3.1661,3.0971,3.2074},
	{2.9223,3.6109,3.6109,3.0123,2.0477,2.8225,2.6882,3.0714,3.4787,3.5695,3.6304,4.038,3.7136,3.3986,3.419,3.3348,3.2784,3.4271,3.1499,3.6221,3.2693,3.6041,3.2023,3.382,3.4965,3.1594,3.1756,3.2063,2.8229,3.1649,3.0123,2.9719,2.7606,3.2234},
	{2.6991,1.665,2.2246,2.7246,2.1122,1.8882,2.1893,2.1269,2.2259,2.3168,2.8094,3.1507,2.5096,2.1794,2.5182,2.1452,1.9975,1.7354,1.7776,1.7632,2.2959,2.2604,2.2513,1.978,2.3479,2.5532,2.5431,2.1947,2.1518,2.1089,2.2391,2.5863,2.547,2.1329},
	{4.714,2.2246,3.6109,3.4177,3.434,3.8341,4.4228,4.3241,3.5657,4.4858,3.4633,4.038,4.5245,4.865,4.5976,4.8752,4.8878,4.6798,6.0403,5.0084,5.0039,4.9904,4.6492,3.8238,3.982,4.4403,3.4439,4.3694,3.4735,3.8581,4.2932,4.1957,4.4188,4.306},
	{2.7681,3.6109,2.5123,4.1109,4.1271,5.2204,4.1352,3.631,3.8841,3.1641,3.3202,2.9882,3.1382,3.191,2.9394,3.0294,3.7892,3.2329,3.1499,3.468,2.5762,2.8309,2.9,3.7184,3.4224,3.3417,3.6109,3.494,3.5605,3.7403,3.8232,4.0134,3.6079,3.3075},
	{2.2717,2.9178,3.6109,4.1109,4.1271,3.6109,3.2189,2.8201,2.4978,2.0655,2.2216,2.2463,2.1266,2.2259,2.1773,2.1452,2.3229,2.0231,1.9972,2.1558,2.1323,2.6878,2.6682,2.6198,2.8034,3.1053,3.1184,3.2708,3.1858,3.635,3.1946,3.0325,3.6079,2.7504},
	{3.2099,3.6109,2.5123,2.3191,3.2108,3.2744,2.9565,2.8201,2.9191,3.3072,3.3202,3.586,3.5129,3.661,3.2113,3.1406,3.907,3.9867,3.9608,4.1611,4.0231,3.6041,3.7329,4.0751,3.4224,3.5648,3.3698,3.0344,2.7003,2.3822,2.7671,2.5863,2.3394,3.0532},
	{5.4072,3.6109,2.9178,4.1109,4.1271,5.2204,4.8283,4.101,4.0177,4.7735,4.642,3.7867,4.8122,5.2704,3.9045,4.0279,4.1947,3.6682,3.5553,3.3344,3.6177,5.3959,4.0896,3.7184,3.8642,3.6518,4.0629,3.8994,4.0125,4.3281,4.1109,5.8051,4.7065,4.6425},
	{3.3277,3.6109,3.6109,4.1109,4.1271,5.2204,4.4228,2.6194,2.5624,2.1109,2.1781,2.6881,2.5786,2.6677,3.2113,3.2658,2.9907,3.3581,3.2677,3.2167,3.1068,3.2558,2.9,3.0766,2.8034,3.1053,3.6109,3.7817,3.5605,3.3726,3.7054,3.8592,3.7257,3.6248},
	{2.574,3.6109,2.9178,2.7246,2.5177,3.4286,2.577,2.6659,2.9678,3.2331,3.9488,3.9045,4.8122,3.3986,3.9045,4.5875,3.5886,3.9867,5.3471,4.7207,3.4635,3.6911,3.956,3.8238,3.7589,4.6634,4.2171,3.6763,3.1252,3.7403,3.5047,3.4072,3.3202,3.0878},
	{3.4613,3.6109,3.6109,3.0123,2.5177,2.33,2.6882,3.0024,3.2555,3.1641,2.6271,2.373,2.5435,2.9191,2.6517,2.8828,2.5525,2.6988,2.6063,2.6105,2.9245,2.6233,2.7033,2.6537,2.4779,2.4662,2.3712,2.29,2.4619,2.5363,2.2917,2.0439,2.7141,2.4602},
	{2.6991,2.9178,3.6109,2.7246,2.8744,2.9178,2.6311,2.7147,2.4978,2.614,2.3394,2.5827,2.6919,2.4372,2.8484,2.9781,2.9419,2.6988,2.7822,2.8881,2.6368,2.5925,2.991,2.93,2.6275,2.3862,2.3979,2.683,2.7804,2.893,2.6068,2.8094,2.6696,2.5672},
	{4.714,3.6109,3.6109,4.1109,4.8203,5.2204,4.1352,5.7104,2.8281,2.8276,3.1379,3.419,3.0204,3.4787,4.1922,4.8752,4.8878,4.9675,4.4308,4.4976,4.4931,3.6911,3.3964,3.0253,3.0657,3.3417,3.6109,4.3694,3.879,4.5512,5.9026,4.1957,3.8592,4.3801},
	{3.6154,3.6109,3.6109,4.1109,4.8203,5.2204,5.5215,3.7645,3.7664,3.7927,2.5019,2.1338,2.2999,2.8725,2.6517,3.2012,3.1532,3.2329,3.0958,2.7748,2.6061,2.2604,1.9924,2.109,2.3003,2.8309,3.4439,3.2708,3.014,2.802,3.13,3.0325,3.4072,2.9377},
	{2.2291,3.6109,3.6109,3.0123,4.8203,5.2204,3.912,4.3241,3.3245,2.8276,2.8502,2.6166,2.3554,2.5624,2.2463,2.3629,2.2729,2.5397,2.4293,2.4961,2.3649,2.4002,2.4801,2.5245,3.1711,2.7539,3.5239,3.581,3.2504,3.4526,3.4177,3.6079,2.9148,3.273}};

int PDB::parse_pdb(string& in, string& top_chain, vector<int>& skip_residues){

	target = in;
	int c = 0;
	char buf[1024];
	std::vector<int>::iterator it;
	string residues[20] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
	map <string, int> aacodes_num; 
	for(int d = 0; d < 20; d++){
		aacodes_num[residues[d]] = d;
	}	

	// Parse PDB file
	FILE *finpdb = fopen (target.c_str(), "r");
	if (finpdb != NULL){
		while (fgets (buf, 1024, finpdb)){
			string line = buf;
			if (strncmp(line.substr(0,4).c_str(),"ATOM",4) == 0){
				all_atoms.push_back(line);
				int tag = 0;
				if(!chains.size()) tag = 1;
				for(unsigned int i = 0; i < chains.size(); i++){
					if (line.substr(21,1) == chains[i]){
						tag = 1;
					}
				}
				if(tag){

					string atom = line.substr(13,4);
					erase_all(atom, " ");
					string res = line.substr(17,3);

					if(strcmp (atom.c_str(),"CA") == 0 && strcmp (res.c_str(),"GLY") == 0){
						double x = atof(line.substr(30,8).c_str());
						double y = atof(line.substr(38,8).c_str());
						double z = atof(line.substr(46,8).c_str());				
						string resnum = line.substr(22,4);
						erase_all(resnum, " ");
						int n = boost::lexical_cast<int>(resnum);
						it = find (skip_residues.begin(), skip_residues.end(), n);	
						if (it == skip_residues.end()){
							Point3d p(aacodes_num[res],n,x,y,z);
							backbone.push_back(p);	
						}else{
							cout << "Skipping residue " << n << endl;
						}
						if(line.substr(21,1) == top_chain){
							pdb_res_nums[n] = c;
						}						
						c++;
					}else if(strcmp (atom.c_str(),"CB") == 0 && strcmp (res.c_str(),"GLY") != 0){
						double x = atof(line.substr(30,8).c_str());
						double y = atof(line.substr(38,8).c_str());
						double z = atof(line.substr(46,8).c_str());		
						string resnum = line.substr(22,4);
						erase_all(resnum, " ");
						int n = boost::lexical_cast<int>(resnum);
						it = find (skip_residues.begin(), skip_residues.end(), n);	
						if (it == skip_residues.end()){
							Point3d p(aacodes_num[res],n,x,y,z);
							backbone.push_back(p);	
						}else{
							cout << "Skipping residue " << n << endl;
						}
						if(line.substr(21,1) == top_chain){
							pdb_res_nums[n] = c;
						}								
						c++;
					}
				}			
			}			
		}
	}else{
		cout << "Couldn't open file " << target << endl;
		return(0);
	}
	fclose(finpdb);
	return(1);
}

void PDB::write_pdb(double energy){

	if(!output.length()){
		output = target;
		output.replace(output.length()-4,10,"_EMBED.pdb");
	}
	ofstream outfile;
	outfile.open (output.c_str());  	

	outfile << "HEADER MEMBRANE_ENERGY " << setprecision(10) << energy << endl;	
	outfile << "HEADER X_ROTATION      " << setprecision(6) << x_rot_opt*(180/M_PI) << endl;
	outfile << "HEADER Y_ROTATION      " << setprecision(6) << y_rot_opt*(180/M_PI) << endl;
	outfile << "HEADER Z_TRANSLATION   " << setprecision(6) << z_trans_opt << endl;

	double max_x = -1000000.0, max_y = -1000000.0, min_x = 1000000.0, min_y = 1000000.0;

	BOOST_FOREACH(string a, all_atoms){

 		int tag = 0;
		if(!chains.size() || output_all_chains){
			tag = 1;
		}else{
			for(unsigned int i = 0; i < chains.size(); i++){
				if (a.substr(21,1) == chains[i]){
					tag = 1;
				}
			}
		}
		if(tag){	
			double x = atof(a.substr(30,8).c_str());
			double y = atof(a.substr(38,8).c_str());
			double z = atof(a.substr(46,8).c_str());	

			Point3d p(0,0,x,y,z);
			transform_atom(p,x_rot_opt,y_rot_opt,z_trans_opt);

			if(p.x_ > max_x) max_x = p.x_;
			if(p.y_ > max_y) max_y = p.y_;	
			if(p.x_ < min_x) min_x = p.x_;
			if(p.y_ < min_y) min_y = p.y_;	

			char buffer[32];
			sprintf(buffer, "%8.3f", p.x_);
			a.replace(30,8,buffer);		
			sprintf(buffer, "%8.3f", p.y_);
			a.replace(38,8,buffer);
			sprintf(buffer, "%8.3f", p.z_);
			a.replace(46,8,buffer);
			outfile << a;
		}
	}

	// Add membrane
	char buffer[1000];
	max_x += 8;
	max_y += 8;
	min_x -= 8;
	min_y -= 8;
	int count = 1000;
	for(int i = 0; i <= (int)(max_x-min_x)/2;i++){
		for(int j = 0; j <= (int)(max_y-min_y)/2;j++){

			double z = 15.000+best_extra_shift;
			double x = min_x + (i * 2);
			double y = min_y + (j * 2);
			snprintf(buffer, 1000, "HETATM %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f",count, " O  ", "DUM", count, x, y, z, 0.0);
			outfile << buffer << endl;
			count++;
			if(count > 9999)count = 1000;
			z = -15.000-best_cyto_shift;
			snprintf(buffer, 1000, "HETATM %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f",count, " N  ", "DUM", count, x, y, z, 0.0);
			outfile << buffer << endl;
			count++;
			if(count > 9999)count = 1000;

			if(polar_flg){
				z = 24.000+best_extra_shift;
				snprintf(buffer, 1000, "HETATM %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f",count, " O  ", "DUM", count, x, y, z, 0.0);
				outfile << buffer << endl;
				count++;
				if(count > 9999)count = 1000;
				z = -24.000 -best_cyto_shift;
				snprintf(buffer, 1000, "HETATM %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f",count, " N  ", "DUM", count, x, y, z, 0.0);
				outfile << buffer << endl;
				count++;
				if(count > 9999)count = 1000;
			}
		}
	}
	outfile.close();
	cout << "Written " << output << endl;

}

int PDB::get_slice_index(double& z){

	// 48  Angstrom thick membrane
	// 32 * 1.5 Angstrom thick slices within the membrane (-24 to 24)
	// 2 slices accounting for cytoplasm/lumen (< -24, > 24)
	// 34 slices

	double slice_start = -24.0;
	double slice_stop = 24.0;

	if(z < slice_start){
		return(0);
	}else if(z >= slice_stop){
		return(33);
	}
	for(int i = 1;i < 33;i++,slice_start+=1.5){
		if(z >= slice_start && z < slice_start+1.5){
			return(i);
		}
	}
	return(0);
}

double PDB::orientate(double xrt, double yrt, double z_trans){

	double energy = 0;
	double xp, yp, zp;
	int top1 = 0;
	int top2 = 0;

	for(unsigned int i = 0; i < backbone.size(); i++){

		double x = backbone[i].x_, y = backbone[i].y_, z = backbone[i].z_;	

		// X-axis rotation
		yp = (y*cos(xrt)) - (z*sin(xrt));
		zp = (y*sin(xrt)) + (z*cos(xrt));
		y = yp;
		z = zp;

		// Y-axis rotation
		zp = (z*cos(yrt)) - (x*sin(yrt));
		xp = (z*sin(yrt)) + (x*cos(yrt));
		z = zp;
		x = xp;

		// Z-axis translation
		z += z_trans;

		if(z > 17.5) top1 = 1;
		if(z < -17.5) top2 = 1;

		energy += mempot[backbone[i].res_][get_slice_index(z)];
	}

	if(force_span && (!top1 || !top2)){
		energy += 10000;
	}
	return(energy);
}

string PDB::get_nterm(double xrt, double yrt, double z_trans){

	string n = "in";
	double yp, zp;
	double x = backbone[0].x_, y = backbone[0].y_, z = backbone[0].z_;	

	// X-axis rotation
	yp = (y*cos(xrt)) - (z*sin(xrt));
	zp = (y*sin(xrt)) + (z*cos(xrt));
	y = yp;
	z = zp;

	// Y-axis rotation
	zp = (z*cos(yrt)) - (x*sin(yrt));
	// Z-axis translation
	z = zp + z_trans;

	if(z > 0){
		n = "out";
	}
	return(n);
}

int PDB::calc_helix_tilt(vector<int>& topology){

	double xp, yp, zp, xq = 0.0, yq = 0.0, zq = 100;
	double top,bottom,angle;

	accumulator_set<double, stats<tag::mean> > x_acc;
	accumulator_set<double, stats<tag::mean> > y_acc;
	accumulator_set<double, stats<tag::mean> > z_acc;
	accumulator_set<double, stats<tag::mean> > z_shift;

	cout << "Calculating helix tilt angles:" << endl << endl;
	cout << "TM Helix\tAngle (degrees)" << endl;

	for(unsigned int i = 0; i < topology.size(); i+= 2){

		cout << topology[i] << "-" << topology[i+1] << "\t\t";

		if(pdb_res_nums.find(topology[i+1]) == pdb_res_nums.end()){
			return(0);
		}
		if(pdb_res_nums.find(topology[i]) == pdb_res_nums.end()){
			return(0);
		}

		z_shift(backbone[pdb_res_nums[topology[i]]].z_);	
		z_shift(backbone[pdb_res_nums[topology[i+1]]].z_);	

		if(backbone[pdb_res_nums[topology[i]]].z_ < backbone[pdb_res_nums[topology[i+1]]].z_){
			xp = backbone[pdb_res_nums[topology[i+1]]].x_ - backbone[pdb_res_nums[topology[i]]].x_;
			yp = backbone[pdb_res_nums[topology[i+1]]].y_ - backbone[pdb_res_nums[topology[i]]].y_;
			zp = backbone[pdb_res_nums[topology[i+1]]].z_ - backbone[pdb_res_nums[topology[i]]].z_;
		}else{
			xp = backbone[pdb_res_nums[topology[i]]].x_ - backbone[pdb_res_nums[topology[i+1]]].x_;
			yp = backbone[pdb_res_nums[topology[i]]].y_ - backbone[pdb_res_nums[topology[i+1]]].y_;
			zp = backbone[pdb_res_nums[topology[i]]].z_ - backbone[pdb_res_nums[topology[i+1]]].z_;
		}

		x_acc(xp);
		y_acc(yp);
		z_acc(zp);
		top = xp*xq + yp*yq + zp*zq;
		bottom = sqrt(xp*xp + yp*yp + zp*zp) * sqrt(xq*xq + yq*yq + zq*zq);
		angle = (180/M_PI)*acos(top/bottom);		
		if(angle > 90) angle = 180-angle;
		cout << setprecision (2) << fixed << angle << endl;

	}
	cout << endl;

	xp = boost::accumulators::mean(x_acc);
	yp = boost::accumulators::mean(y_acc);
	zp = boost::accumulators::mean(z_acc);
	double mean_tm_z = boost::accumulators::mean(z_shift);
	top = xp*xq + yp*yq + zp*zq;
	bottom = sqrt(xp*xp + yp*yp + zp*zp) * sqrt(xq*xq + yq*yq + zq*zq);
	angle = (180/M_PI)*acos(top/bottom);		
	if(angle > 90) angle = 180-angle;

	cout << "Average tilt:\t\t" << setprecision (2) << fixed << angle << endl;
	cout << "Average TM Z-coord:\t" << setprecision (2) << fixed << mean_tm_z << endl;

	return(1);
}


int PDB::pre_position(string& n_term, vector<int>& topology){

	double xp,yp,zp,xq = 0.0,yq = 0.0,zq = 100;

	accumulator_set<double, stats<tag::mean> > x_acc;
	accumulator_set<double, stats<tag::mean> > y_acc;
	accumulator_set<double, stats<tag::mean> > z_acc;
	accumulator_set<double, stats<tag::mean> > z_shift;

	int toggle = 0;
	for(unsigned int i = 0; i < topology.size(); i+= 2){

		if(pdb_res_nums.find(topology[i+1]) == pdb_res_nums.end()){
			return(0);
		}
		if(pdb_res_nums.find(topology[i]) == pdb_res_nums.end()){
			return(0);
		}

		if(toggle){

			xp = backbone[pdb_res_nums[topology[i+1]]].x_ - backbone[pdb_res_nums[topology[i]]].x_;
			yp = backbone[pdb_res_nums[topology[i+1]]].y_ - backbone[pdb_res_nums[topology[i]]].y_;
			zp = backbone[pdb_res_nums[topology[i+1]]].z_ - backbone[pdb_res_nums[topology[i]]].z_;
			x_acc(xp);
			y_acc(yp);
			z_acc(zp);
			toggle = 0;

		}else{

			xp = backbone[pdb_res_nums[topology[i]]].x_ - backbone[pdb_res_nums[topology[i+1]]].x_;
			yp = backbone[pdb_res_nums[topology[i]]].y_ - backbone[pdb_res_nums[topology[i+1]]].y_;
			zp = backbone[pdb_res_nums[topology[i]]].z_ - backbone[pdb_res_nums[topology[i+1]]].z_;
			x_acc(xp);
			y_acc(yp);
			z_acc(zp);
			toggle = 1;
		}	

	}

	double xp_a = boost::accumulators::mean(x_acc);
	double yp_a = boost::accumulators::mean(y_acc);
	double zp_a = boost::accumulators::mean(z_acc);

	double top = xp_a*xq + yp_a*yq + zp_a*zq;
	double bottom = sqrt(xp_a*xp_a + yp_a*yp_a + zp_a*zp_a) * sqrt(xq*xq + yq*yq + zq*zq);
	y_rot_pre = acos(top/bottom);
	if(xp_a > 0) y_rot_pre *= -1;
	z_rot_pre = -atan(yp_a/xp_a);

	// Rotate backbone residues
	BOOST_FOREACH(Point3d & p, backbone){
		double xi = (p.x_*cos(z_rot_pre) - p.y_*sin(z_rot_pre));
		double yi = (p.x_*sin(z_rot_pre) + p.y_*cos(z_rot_pre)); 
		p.x_ = xi;
		p.y_ = yi;

		double zi = (p.z_*cos(y_rot_pre)) - (p.x_*sin(y_rot_pre));
		xi = (p.z_*sin(y_rot_pre)) + (p.x_*cos(y_rot_pre));
		p.x_ = xi;
		p.z_ = zi;
	}

	// Calculate mean Z-coord for TM residues
	for(unsigned int i = 0; i < topology.size(); i+= 2){
		if(!i){
			if(backbone[pdb_res_nums[topology[i]]].z_ < backbone[pdb_res_nums[topology[i+1]]].z_){
				if(n_term == "out") pre_flip = true;
			}else{
				if(n_term == "in") pre_flip = true;
			}
		}
		z_shift(backbone[pdb_res_nums[topology[i]]].z_);	
		z_shift(backbone[pdb_res_nums[topology[i+1]]].z_);	
	}

	// Translate Z-coords so helices are centred on Z=0
	z_trans_pre = boost::accumulators::mean(z_shift);
	BOOST_FOREACH(Point3d & p, backbone){
		p.z_ -= z_trans_pre;
		// Flip if N-terminal is wrong
		if(pre_flip){
			double yp = (p.y_*cos(M_PI)) - (p.z_*sin(M_PI));
			double zp = (p.y_*sin(M_PI)) + (p.z_*cos(M_PI));
			p.y_ = yp;
			p.z_ = zp;			
		}
	}	
	pre_pos = true;
	return(1);
}

void PDB::transform_atom(Point3d& p, double xrt, double yrt, double z_trans){
	
	double xp, yp, zp;
	double x = p.x_, y = p.y_, z = p.z_;

	x -= max_x;
	y -= max_y;
	z -= max_z;

	if(pre_pos){
		double xi = (x*cos(z_rot_pre) - y*sin(z_rot_pre));
		double yi = (x*sin(z_rot_pre) + y*cos(z_rot_pre)); 
		x = xi;
		y = yi;

		double zi = (z*cos(y_rot_pre)) - (x*sin(y_rot_pre));
		xi = (z*sin(y_rot_pre)) + (x*cos(y_rot_pre));
		x = xi;
		z = zi;

		z -= z_trans_pre;

		// Flip if N-terminal is wrong
		if(pre_flip){
			double yp = (y*cos(M_PI)) - (z*sin(M_PI));
			double zp = (y*sin(M_PI)) + (z*cos(M_PI));
			y = yp;
			z = zp;			
		}
	}

	// X-axis rotation	
	yp = (y*cos(xrt)) - (z*sin(xrt));
	zp = (y*sin(xrt)) + (z*cos(xrt));
	y = yp;
	z = zp;

	// Y-axis rotation
	zp = (z*cos(yrt)) - (x*sin(yrt));
	xp = (z*sin(yrt)) + (x*cos(yrt));
	z = zp;
	x = xp;
	
	z += z_trans;

	if(flip){
		yp = (y*cos(M_PI)) - (z*sin(M_PI));
		zp = (y*sin(M_PI)) + (z*cos(M_PI));
		y = yp;
		z = zp;
	}

	p.x_ = x;
	p.y_ = y;
	p.z_ = z;

}

void PDB::calc_thickness(double x_rot, double y_rot, double z_trans){

	double xrt = x_rot;
	double yrt = y_rot;	
	double xp, yp, zp;
	double lowest_e = 100000;
	double shift_increment = 0.25;
	double aa_average[20];     // average in membrane core
	double aa_average_out[20]; // average outside membrane
	double head_top = 20.0;
	double head_bot = 10.0;
	double shift_start = -7.5;	
	double shift_stop = 10.0;	
	double shift_start2 = -7.5;	
	double shift_stop2 = 10.0;	

	if(beta){
		shift_start = -10.0;	
		shift_stop = -4.0;	
		shift_start2 = -10.0;	
		shift_stop2 = -2.0;	
	}

	//Average for core region
	for(int i = 0; i < 20; i++){
		double ave = 0;
		for(int j = get_slice_index(shift_start); j <= get_slice_index(shift_stop); j++){
			ave += mempot[i][j];
		}
		aa_average[i] = ave/(1 + get_slice_index(shift_stop) - get_slice_index(shift_start));
	}

	//Average outside membrane
	for(int i = 0; i < 20; i++){
		double ave = 0;
		int count = 0;
		for(int j = 0; j < get_slice_index(shift_start); j++){
			ave += mempot[i][j];
			count++;
		}
		for(int j = get_slice_index(shift_stop); j <= 33; j++){
			ave += mempot[i][j];
			count++;
		}
		aa_average_out[i] = ave/count;
	}

	// Get max and min Z coords post transformation
	double max_zt = -10000;
	double min_zt =  10000;
	BOOST_FOREACH(Point3d p, backbone){
		double x = p.x_;
		double y = p.y_;
		double z = p.z_;
		yp = (y*cos(xrt)) - (z*sin(xrt));
		zp = (y*sin(xrt)) + (z*cos(xrt));
		y = yp;
		z = zp;
		zp = (z*cos(yrt)) - (x*sin(yrt));
		xp = (z*sin(yrt)) + (x*cos(yrt));
		z = zp;
		z += z_trans;
		if(z > max_zt){
			max_zt = z;
		}
		if(z < min_zt){
			min_zt = z;
		}
	}

	for(double shift_e = shift_start; shift_e < shift_stop; shift_e += shift_increment){
		double energy_e = 0;
		BOOST_FOREACH(Point3d p, backbone){

			double z = p.z_;
			if(xrt||yrt||z_trans){
				double x = p.x_;
				double y = p.y_;
				// X-axis rotation
				yp = (y*cos(xrt)) - (z*sin(xrt));
				zp = (y*sin(xrt)) + (z*cos(xrt));
				y = yp;
				z = zp;

				// Y-axis rotation
				zp = (z*cos(yrt)) - (x*sin(yrt));
				xp = (z*sin(yrt)) + (x*cos(yrt));
				z = zp;
				x = xp;

				// Z-axis translation
				z += z_trans;
			}
			z -= shift_e;

			if(z >= head_bot && z <= head_top){	
				energy_e += mempot[p.res_][get_slice_index(z)];
			}else if(z < head_bot){
				energy_e += aa_average[p.res_];
			}else{
				energy_e += aa_average_out[p.res_];
			}
		}
		
		for(double shift_c = shift_start2; shift_c < shift_stop2; shift_c += shift_increment){

			double energy_c = 0;
			BOOST_FOREACH(Point3d q, backbone){

				double z = q.z_;
				if(xrt||yrt||z_trans){
					double x = q.x_;
					double y = q.y_;
					// X-axis rotation
					yp = (y*cos(xrt)) - (z*sin(xrt));
					zp = (y*sin(xrt)) + (z*cos(xrt));
					y = yp;
					z = zp;

					// Y-axis rotation
					zp = (z*cos(yrt)) - (x*sin(yrt));
					xp = (z*sin(yrt)) + (x*cos(yrt));
					z = zp;
					x = xp;

					// Z-axis translation
					z += z_trans;
				}

				double zx = z + shift_c;
				if(zx <= -head_bot && zx >= -head_top){	
					energy_c += mempot[q.res_][get_slice_index(zx)];
				}else if (zx > -head_top){
					energy_c += aa_average[q.res_];
				}else{
					energy_c += aa_average_out[q.res_];
				}
			}

			
			double total = energy_e+energy_c;
			
			if((total < lowest_e)&&(15+best_extra_shift < max_zt && (-15-best_cyto_shift) > min_zt)){
				/*	
				cout << "thickness = " << (30.0+best_extra_shift+best_cyto_shift) << " - cyto: " << best_cyto_shift << " extra: " << best_extra_shift;
				cout << endl << 15+best_extra_shift << " 0===== =====0 " << -15-best_cyto_shift << endl;				
				cout << max_zt << " <========> " << min_zt << endl;
				cout << "Energy: " << total << endl << endl;
				*/
				lowest_e = total;
				best_extra_shift = shift_e;
				best_cyto_shift = shift_c;

			}
		}
	}
	printf("Hydrophobic thickness:\t%f\n\n",30.0+best_extra_shift+best_cyto_shift);


}

void PDB::origin_shift(){

    BOOST_FOREACH(Point3d p, backbone){

		if(p.x_ > max_x) max_x = p.x_;
		if(p.y_ > max_y) max_y = p.y_;
    	if(p.z_ > max_z) max_z = p.z_;

		BOOST_FOREACH(Point3d q, backbone){
			double dx = p.x_ - q.x_;
			double dy = p.y_ - q.y_;	
			double dz = p.z_ - q.z_;
			double dist = dx*dx + dy*dy + dz*dz;				
			if(dist > max_c_dist) max_c_dist = dist;
		}	
	}
	max_c_dist = sqrt(max_c_dist);

	BOOST_FOREACH(Point3d & p, backbone){
		p.x_ -= max_x;	
		p.y_ -= max_y;	
		p.z_ -= max_z;	
	}	
}

void PDB::set_beta(bool t){

	cout << "Using beta-barrel potential...\n";
	beta = t;
	for(int i = 0; i < 20; i++){
		for(int j = 0; j < 34; j++){
			mempot[i][j] = betamempot[i][j];
		}		
	}
}

int PDB::parse_potential(string& mempotfile){

	int c = 0;
	char buf[1024];
	if(mempotfile.length()){
		FILE *fin = fopen (mempotfile.c_str(), "r");
		if (fin != NULL){
			while (fgets (buf, 1024, fin)){
				string line = buf;
				if (sscanf(buf, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&mempot[c][0],&mempot[c][1],&mempot[c][2],&mempot[c][3],&mempot[c][4],&mempot[c][5],&mempot[c][6],&mempot[c][7],&mempot[c][8],&mempot[c][9],&mempot[c][10],&mempot[c][11],&mempot[c][12],&mempot[c][13],&mempot[c][14],&mempot[c][15],&mempot[c][16],&mempot[c][17],&mempot[c][18],&mempot[c][19],&mempot[c][20],&mempot[c][21],&mempot[c][22],&mempot[c][23],&mempot[c][24],&mempot[c][25],&mempot[c][26],&mempot[c][27],&mempot[c][28],&mempot[c][29],&mempot[c][30],&mempot[c][31],&mempot[c][32],&mempot[c][33]) == 34){
					c++;	
				}
			}
		}else{
			cout << "Couldn't open file " << mempotfile << endl;
			return(0);
		}
		if(c != 20){
			cout << "Failed to parse " << mempotfile << endl;
			return(0);
		}else{
			cout << "Using potential file " << mempotfile << endl;
		}

		if(fin != NULL){
			fclose(fin);	
		}
	}
	return(1);

}