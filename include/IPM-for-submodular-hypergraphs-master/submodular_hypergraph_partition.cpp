//This version only supports unweighted + card-based 
#include <math.h>
#include <stdio.h>
#include <string>
#include <random>
#include "mex.h"
#include "matrix.h"

double mean(double * data, int len){
	double sum = 0;
	for(int i = 0; i < len; i++){
		sum += data[i];
	}
	return sum/len;
}
void swap_double(double * a, int p, int q){
	double temp = a[p];
	a[p] = a[q];
	a[q] = temp;
}
void swap_int(int * a, int p, int q){
	int temp = a[p];
	a[p] = a[q];
	a[q] = temp;
}
void MYsort(double * a, double * b, int len, std::string type){
	// sort a according to rule type, permutate b according to a 
	if(len <= 1) return;
 	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, len-1);
	int pos = distribution(generator), p = 0;
	double pivot = a[pos];
    swap_double(a, pos, len-1);
    swap_double(b, pos, len-1);
	if(type == "descend"){
		for(int i = 0; i < len-1; i++){
			if(a[i] > pivot){
				swap_double(a, i, p);
				swap_double(b, i, p);
			    p++;
			}
		}
		swap_double(a, p, len-1);
		swap_double(b, p, len-1);
	}
	MYsort(a, b, p, type);
	MYsort(a+p+1, b+p+1, len-p-1, type);
}
void INVindex(double * invindex, double * index, int len){
	for(int i = 0; i < len; i++){
		invindex[(int)index[i]] = i;
	}
}

void remove_median (double * x, double * mu, int len){
    double * mu_index = new double [len];
    double * xsort = new double [len];
    int i, pos;
    for(i = 0; i < len; i++){
        mu_index[i] = mu[i];
        xsort[i] = x[i];
    }
    MYsort(xsort, mu_index, len, "descend");
    for(i = 1; i < len; i++){
        mu_index[i] += mu_index[i-1];
    }
    for(i = 0; i < len; i++){
        if(mu_index[i] >= mu_index[len-1]/2){
            pos = i;
            break;
        }
    }
    for(i = 0; i < len; i++){
        x[i] -= xsort[pos];
    }
}

void concave_card_subroutine_Noweight(double * x, double lambdamin, double lambdamax, double * a, double * para, int i, int len){
	int j;
	double minval = INFINITY, tempsum = 0, lambda, val;
	int minindex;
	if(i%2 == 0){
		lambda = mean(a, len) - mean(para, len);
		for(j = 0; j < len; j++){
			val = lambda - a[j];
			tempsum = tempsum + val + para[j];
			if(tempsum < minval){
				minval = tempsum;
				minindex = j;
			}
		}
		if(minval > 0 || minindex == len - 1){
			for(j = 0; j < len; j++){
				x[j] = lambda;
			}
			return;
		}
	}else{
		lambda = (lambdamin + lambdamax)/2;
		for(j = 0; j < len; j++){
			val = lambda - a[j];
			tempsum = tempsum + val + para[j];
			if(tempsum < minval){
				minval = tempsum;
				minindex = j;
			}
		}
		if(minval > 0){
			concave_card_subroutine_Noweight(x, lambdamin, lambda, a, para, i+1, len);
			return;
		}
		if(minindex == len-1){
			concave_card_subroutine_Noweight(x, lambda, lambdamax, a, para, i+1, len);
			return;
		}		
	}
	concave_card_subroutine_Noweight(x, lambda, lambdamax, a, para, i+1, minindex+1);
	concave_card_subroutine_Noweight(x + minindex+1, lambdamin, lambda, a + minindex + 1, para + minindex + 1, i+1, len - minindex - 1);
}

void projection_concave_card_Noweight(double * y, double * a, double * para, int len){
	double * index = new double[len];
    double * asort = new double[len];
	for(int i = 0; i < len; i++){
		index[i] = i;
        asort[i] = a[i];
	}
	double * invindex = new double[len];
	double lambdamin;
	double lambdamax;
	MYsort(asort, index, len, "descend");
	INVindex(invindex, index, len);
	lambdamin = - para[0] + asort[len-1];
	lambdamax = - para[len-1] + asort[0];
	double * x = new double[len]();
	concave_card_subroutine_Noweight(x, lambdamin, lambdamax, asort, para, 0, len);
    for(int i = 0; i < len; i++){
		y[i] = asort[(int)invindex[i]] - x[(int)invindex[i]];
	}
	delete[] index;
	delete[] invindex;
	delete[] x;
    delete[] asort;
}

double eval_Q1(double * x, int ** incidence_list, int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size, int N, int R){
    double Q1 = 0;
    double * tempx;
    double * index;
    int i, j;
    for(i = 0; i < R; i++){
        //TODO: check submodular type
        tempx = new double [incidence_list_size[i]];
        index = new double [incidence_list_size[i]]();
        for(j = 0; j < incidence_list_size[i]; j++){
            tempx[j] = x[incidence_list[i][j]];
        }
        MYsort(tempx, index, incidence_list_size[i], "descend");
        for(j = 0; j < incidence_list_size[i]; j++){
            Q1 += tempx[j]*parameter_list[i][j];
        }
        delete[] tempx;
        delete[] index;
    }
    return Q1;
}

double eval_balanced(double * x, int ** incidence_list, int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size, int N, int R,
        double * mu){
    // typically require remove_median of x
    double Q1, norm_mu_x = 0;
    Q1 = eval_Q1(x, incidence_list, incidence_list_size,
                        parameter_list, parameter_list_size,
                        submodular_type, submodular_type_size, N, R);
    for(int i = 0;i < N; i++){
        norm_mu_x += mu[i]*fabs(x[i]);
    }
    return Q1/norm_mu_x;
}

void optthreshold(double * labels, double * incidence_partition,
        double * NCut, double * x, int ** incidence_list,
        int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size,
        double * mu, int N, int R){
    double * index = new double [N];
    double * xsort = new double [N];
    double * xtemp = new double [N]();
    int i, j, optpos;
    double Q1, tempvol = 0, sum_mu = 0, tempval;
    *NCut = INFINITY;
    for(i = 0; i < N; i++){
        index[i] = i;
        xsort[i] = x[i];
        sum_mu += mu[i];
    }
    MYsort(xsort, index, N, "descend");
    for(i = 0; i < N-1; i++){
        xtemp[(int)index[i]] = 1;
        Q1 = eval_Q1(xtemp, incidence_list, incidence_list_size,
                parameter_list, parameter_list_size,
                submodular_type, submodular_type_size, N, R);
        tempvol += mu[(int)index[i]];
        tempval = Q1/fmin(tempvol, sum_mu - tempvol);
        if(tempval < *NCut){
            *NCut = tempval;
            optpos = i;
        }
    }
    for(i = 0; i < N; i++){
        if(i <= optpos){
            labels[(int)index[i]] = 1;
        }else{
            labels[(int)index[i]] = 0;
        }
    }
    for(i = 0; i < R; i++){
        incidence_partition[i] = labels[incidence_list[i][0]];
        for(j = 1; j < incidence_list_size[i]; j++){
            if(fabs(incidence_partition[i] - labels[incidence_list[i][j]]) > 1e-6){
                incidence_partition[i] = 2;
                break;
            }
        }
    }
    delete[] xsort;
    delete[] xtemp;
    delete[] index;
}

void derivative_mu_norm(double * x, double * g, double * mu, int N){
    double mu_x_pos=0;
    double mu_x_neg=0;
    double mu_x_zeros=0;
    double muzero;
    double epsilon = 1e-10;
    int i;
    for(i = 0; i < N; i++){
        if(x[i] > epsilon){
            mu_x_pos += mu[i];
        }
        if(x[i] < - epsilon){
            mu_x_neg += mu[i];
        }
        if(x[i] < epsilon && x[i] > - epsilon){
            mu_x_zeros += mu[i];
        }
    }
    muzero = -(mu_x_pos - mu_x_neg)/mu_x_zeros;
    for(i = 0; i < N; i++){
        if(x[i] > epsilon){
            g[i] = mu[i];
        }
        if(x[i] < - epsilon){
            g[i] = -mu[i];
        }
        if(x[i] < epsilon && x[i] > - epsilon){
            g[i] = muzero * mu[i];
        }
    }    
};


void main_func(double * labels, double * incidence_partition,
        double * NCut, double *gap,
        int ** incidence_list, int * incidence_list_size,
        double ** parameter_list, int * parameter_list_size,
        char ** submodular_type, int * submodular_type_size,
        double * mu, int N, int R,
        double dec_outloop, double err_inloop, double * warmstart){
    double * x = new double [N]();
    double * x_final = new double [N]();
    double * g = new double [N]();
    unsigned int T = 10000*R;
    int record_dis = (T+99)/100;
    double ** y = new double * [R];
    double ** newy = new double * [R];
    double ** a = new double * [R];
    double * y_sum = new double [N]();
    int i, k, t;
    double tempeta, eta;
 	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, R-1);
    int picked, counter=0;
    double dv, normx, Q1, pv, xdotg, norm_mu_x;
    for(i = 0; i < N; i++){
        x[i] = warmstart[i];
    }
    for(i = 0; i < R; i++){
         //mexPrintf("%d\n", incidence_list_size[i]);
        y[i] = new double [incidence_list_size[i]]();
        newy[i] = new double [incidence_list_size[i]]();
        a[i] = new double [incidence_list_size[i]]();
        projection_concave_card_Noweight(y[i], a[i], parameter_list[i], incidence_list_size[i]);
        for(k = 0; k < incidence_list_size[i]; k++){
            y_sum[incidence_list[i][k]] += y[i][k];
        }
    }

    remove_median(x, mu, N);
    tempeta = eval_balanced(x, incidence_list, incidence_list_size,
            parameter_list, parameter_list_size,
            submodular_type, submodular_type_size, N, R, mu);
    
    while(1){
        eta = tempeta;
        derivative_mu_norm(x, g, mu, N);
        for(t = 0; t < T; t++){
            picked = distribution(generator);
            for(i = 0; i < incidence_list_size[picked]; i++){
                a[picked][i] = y[picked][i] - y_sum[incidence_list[picked][i]] + eta*g[incidence_list[picked][i]];
            }
            projection_concave_card_Noweight(newy[picked], a[picked], parameter_list[picked], incidence_list_size[picked]);
            
            for(i = 0; i < incidence_list_size[picked]; i++){
                y_sum[incidence_list[picked][i]] += newy[picked][i] - y[picked][i];
                y[picked][i] = newy[picked][i];
            }
            if((t+1)%record_dis == 0){
                normx = 0;
                for(i = 0; i < N; i++){
                    x[i] = eta*g[i] - y_sum[i];
                    normx += x[i]*x[i];
                }
                normx = sqrt(normx);
                dv = -normx;
                for(i = 0; i < N; i++){
                    x[i] = x[i]/normx;
                }
                Q1 = eval_Q1(x, incidence_list, incidence_list_size,
                        parameter_list, parameter_list_size,
                        submodular_type, submodular_type_size, N, R);
                xdotg = 0;
                for(i = 0; i < N; i++){
                    xdotg += x[i]*g[i];
                }
                pv = Q1 - eta*xdotg;
                remove_median(x, mu, N);
                norm_mu_x = 0;
                for(i = 0; i < N; i++){
                    norm_mu_x += mu[i]*fabs(x[i]);
                }
                tempeta = Q1/norm_mu_x;
                *gap = pv - dv;
                mexPrintf("%d, gap: %f, tempeta: %f\n", t, *gap, tempeta);
                if(tempeta < (1-dec_outloop-1e-5)*eta || *gap < err_inloop){
                    if(tempeta < eta){
                        for(i = 0; i < N; i++){
                            x_final[i] = x[i];
                        }
                    }
                    break;
                }
            }
        }
        if(tempeta > (1-dec_outloop)*eta){
            break;
        }
    }
    optthreshold(labels, incidence_partition, NCut, x_final,
            incidence_list, incidence_list_size,
            parameter_list, parameter_list_size,
            submodular_type, submodular_type_size, mu, N, R);
    delete[] x;
    delete[] x_final;
    delete[] g;
    delete[] y;
    delete[] newy;
    delete[] a;
    delete[] y_sum;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	if (nlhs != 4 || nrhs != 9) {
    mexWarnMsgTxt("Check Parameters");
    return;
    }
    // read incidence_list
    int R = *(mxGetPr(prhs[5]));
    const mxArray * incidence_list_org = prhs[0];
    mxArray * incidence_Element;
    //const int * R_org = mxGetDimensions(incidence_list_org);
    int * * incidence_list = new int * [R];
    int * incidence_list_size = new int [R];
    double * templist;
    int j, k;
    mexPrintf("%d\n", 1);
    for(j = 0; j < R;j++){
       incidence_Element = mxGetCell(incidence_list_org, j);
       incidence_list_size[j] = (int)mxGetN(incidence_Element);
       incidence_list[j] = new int[incidence_list_size[j]];
       templist = mxGetPr(incidence_Element);
       for(k = 0; k < incidence_list_size[j]; k++){
           incidence_list[j][k] = (int)templist[k]-1;
           //mexPrintf("%d, %d, %d\n", j, k, incidence_list[j][k]);
       }
    }
    // read parameter_list
    const mxArray * parameter_list_org = prhs[1];
    mxArray * parameter_Element;
    double * * parameter_list = new double * [R];
    int * parameter_list_size = new int [R];
    for(j = 0; j < R;j++){
       parameter_Element = mxGetCell(parameter_list_org, j);
       parameter_list_size[j] = (int)mxGetN(parameter_Element);
       parameter_list[j] = new double[parameter_list_size[j]];
       templist = mxGetPr(parameter_Element);
       for(k = 0; k < parameter_list_size[j]; k++){
           parameter_list[j][k] = templist[k];
           //mexPrintf("%d, %d, %f\n", j, k, parameter_list[j][k]);
       }
    }
    // read submodular_type
    const mxArray * submodular_type_org = prhs[2];
    mxArray * submodular_type_Element;
    char ** submodular_type = new char * [R];
    int * submodular_type_size = new int [R];
    char * temp_type;
    for(j = 0; j < R;j++){
       submodular_type_Element = mxGetCell(submodular_type_org, j);
       temp_type = (char *)mxGetPr(submodular_type_Element);
       submodular_type_size[j] = (int)mxGetN(submodular_type_Element);
       submodular_type[j] = new char[submodular_type_size[j]];
       for(k=0; k < submodular_type_size[j]; k++){
           submodular_type[j][k] = (char)temp_type[2*k];
           //mexPrintf("%d, %d, %c\n", j, k, temp_type[2*k]);
       }
    }
    
    double * mu = mxGetPr(prhs[3]);
    int N = *(mxGetPr(prhs[4]));
    double dec_outloop = *(mxGetPr(prhs[6]));
    double err_inloop = *(mxGetPr(prhs[7]));
    double * warmstart = mxGetPr(prhs[8]);
    
    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    double * labels = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(R, 1, mxREAL);
    double * incidence_partition = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double * NCut = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double * gap = mxGetPr(plhs[3]);
    
    main_func(labels, incidence_partition, NCut, gap,
            incidence_list, incidence_list_size,
            parameter_list, parameter_list_size,
            submodular_type, submodular_type_size, 
            mu, N, R, dec_outloop,
            err_inloop, warmstart);
    delete[] incidence_list;
    delete[] incidence_list_size;
    delete[] parameter_list_size;
    delete[] parameter_list;
    delete[] submodular_type;
    delete[] submodular_type_size;
}







