#include <stdio.h>
#include <algorithm>
#include <cstring>
#include <map>
#include <iostream>
#include <set>
#include <cmath>
#include <vector>
#include <time.h>
#include <pthread.h>
#include <atomic>

using namespace std;

#define R_A 0.001388889/5
#define R_B 0.001388889/5
#define DIM 2
#define THREADNUM 4096

const double h = M_PI/180.0;
const double DIS = 9.0*(R_A*R_A+R_B*R_B);
const double arcDis = sqrt(DIS)*M_PI/180.0;
const double eDis = (2*(1 - cos(sqrt(DIS)*h)));

typedef struct node{
    double x, y;
}node;

typedef struct BuildArg{
	unsigned int l, r, k;
    BuildArg(unsigned int a, unsigned int b, unsigned int c):l(a),r(b),k(c){}
}BuildArg;

typedef struct ThreadArg{
	unsigned int l, r, x, xMax;
}ThreadArg;

atomic<unsigned long long> prume = 0;

node *s, *s2;
unsigned int *d, *lc, *rc;
double *L, *R, *D, *U;
unsigned int n;
int queryThreadNum = 4096, curThreadNum = 512;



inline double dist1(unsigned int a, unsigned int b) { // �� 
  return (s[a].x - s[b].x) * (s[a].x - s[b].x) +
         (s[a].y - s[b].y) * (s[a].y - s[b].y);
}

inline double dist2(unsigned int a, unsigned int b){ // 
  return (s[a].x - s2[b].x) * (s[a].x - s2[b].x) * cos((s[a].y + s2[b].y)/360*M_PI) +
         (s[a].y - s2[b].y) * (s[a].y - s2[b].y);
}

inline double dist3(unsigned int a, unsigned int b){ // , second is s2
	//return (acos(cos(s[a].y)*cos(s[b].y)*cos(s[a].x-s[b].x) + sin(s[a].y)*sin(s[b].y)));
	return 2*(asin(sqrt(sin((s[a].y-s2[b].y)/2)*sin((s[a].y-s2[b].y)/2) + cos(s[a].y)*cos(s2[b].y)*sin((s[a].x-s2[b].x)/2)*sin((s[a].x-s2[b].x)/2))));
}

double dist4(unsigned int a, unsigned int b){
	double x1, x2, y1, y2, z1, z2;

	x1 = cos(s[a].x)*sin((M_PI_2-s[a].y));
	y1 = sin(s[a].x)*sin((M_PI_2-s[a].y));
	z1 = cos((M_PI_2-s[a].y));

	x2 = cos(s2[b].x)*sin((M_PI_2-s2[b].y));
	y2 = sin(s2[b].x)*sin((M_PI_2-s2[b].y));
	z2 = cos((M_PI_2-s2[b].y));
	
	return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);

}

bool cmp1(node a, node b) { return a.x < b.x;}

bool cmp2(node a, node b) { return a.y < b.y;}

void maintain(int x){
    L[x] = R[x] = s[x].x;
    D[x] = U[x] = s[x].y;
    if (lc[x]){
        L[x] = min(L[x], L[lc[x]]), R[x] = max(R[x], R[lc[x]]);
        D[x] = min(D[x], D[lc[x]]), U[x] = max(U[x], U[lc[x]]);
    }
    if (rc[x]){
        L[x] = min(L[x], L[rc[x]]), R[x] = max(R[x], R[rc[x]]);
        D[x] = min(D[x], D[rc[x]]), U[x] = max(U[x], U[rc[x]]);
    }
}


// 维度轮换

unsigned int build(unsigned int l, unsigned int r, unsigned int sdim){
	
    if (l > r) return 0;
    if (l == r){
        maintain(l);
        return l;
    }

    unsigned int mid = (l + r)/2;

    if (sdim == 0){
        nth_element(s+l, s+mid, s+r+1, cmp1);
    }else{
        nth_element(s+l, s+mid, s+r+1, cmp2);
    }

    lc[mid] = build(l, mid-1, (sdim+1)%2), rc[mid] = build(mid+1, r, (sdim+1)%2);
    maintain(mid);

    return mid;
}


void* threadBuild(void* arg){
	unsigned int l, r, k;
	l = ((BuildArg*)arg)->l;
	r = ((BuildArg*)arg)->r;
	k = ((BuildArg*)arg)->k;
	
	unsigned int* result = (unsigned int*)malloc(sizeof(unsigned int));
	*result = 0;
	
    if (l > r) {
		pthread_exit(result);
	}
    if (l == r){
        maintain(l);
		*result = l;
		pthread_exit(result);
    }

    unsigned int mid = (l+r) >> 1;

    if (k == 0){
        nth_element(s+l, s+mid, s+r+1, cmp1);
    }else if (k == 1){
        nth_element(s+l, s+mid, s+r+1, cmp2);
    }

	if (r-l >= n/curThreadNum){
		int res;
		void* threadRet = NULL;
		
		pthread_t threadL, threadR;
		BuildArg threadRArgs = {mid+1, r, (k+1)%2};
		BuildArg threadLArgs = {l, mid-1, (k+1)%2};
		
		res = pthread_create(&threadL, NULL, threadBuild, (void*)(&threadLArgs));
		if (res != 0){
			printf("Create thread failed");
			exit(res);
		}
		res = pthread_create(&threadR, NULL, threadBuild, (void*)(&threadRArgs));
		if (res != 0){
			printf("Create thread failed");
			exit(res);
		}

		// wait to join
		res = pthread_join(threadL, &threadRet);
		if (res){
			printf("thread not joined.\n");
		}
		lc[mid] = *((unsigned int*)(threadRet));
		free(threadRet);
		
		res = pthread_join(threadR, &threadRet);
		if (res){
			printf("thread not joined.\n");
		}
		rc[mid] = *((unsigned int*)(threadRet));
		free(threadRet);

		maintain(mid);

	}else{
    	lc[mid] = build(l, mid-1, (k+1)%2), rc[mid] = build(mid+1, r, (k+1)%2);
	    maintain(mid);

	}

	*result = mid;
	pthread_exit(result);
}






unsigned int build(unsigned int l, unsigned int r){
    if (l > r) return 0;
    if (l == r){
        maintain(l);
        return l;
    }

    unsigned int mid = (l+r) >> 1;
    double avx = 0, avy = 0, vax = 0, vay = 0;
    for (unsigned int i=l; i<=r; i++){
        avx += s[i].x, avy += s[i].y;
    }
    avx /= (double)(r - l + 1);
    avy /= (double)(r - l + 1);

    for (int i=l; i<=r; i++){
        vax += (s[i].x - avx) * (s[i].x - avx);
        vay += (s[i].y - avy) * (s[i].y - avy);
    }

    if (vax >= vay){
        d[mid] = 1;
        nth_element(s+l, s+mid, s+r+1, cmp1);
    }else{
        d[mid] = 2;
        nth_element(s+l, s+mid, s+r+1, cmp2);
    }

    lc[mid] = build(l, mid-1), rc[mid] = build(mid+1, r);
    maintain(mid);

    return mid;
}

double f(unsigned int b, unsigned int a){ // a is s2
    double ret = 0;
    double ra, dec, x, y, alpha, cosTheta, tmp;
    

    if (!(s2[a].x >= L[b] && s2[a].x <= R[b])){
		ra = s2[a].x;
		dec = M_PI_2-s2[a].y;
		x = cos(ra)*sin(dec), y = sin(ra)*sin(dec);
		alpha = L[b];
		cosTheta = (cos(alpha+M_PI_2)*x + sin(alpha+M_PI_2)*y);
		ret = asin(abs(cosTheta));
		alpha = R[b];
		cosTheta = (cos(alpha+M_PI_2)*x + sin(alpha+M_PI_2)*y);
		ret = min(ret, asin(abs(cosTheta)));
	}
	
	if (!(s2[a].y >= D[b] && s2[a].y <= U[b])){
		ra = s2[a].x, dec = s2[a].y;
		tmp = abs(dec - D[b]);
		if (ret == 0) ret = tmp;
		else ret = min(ret, tmp);
		tmp = abs(U[b] - dec);
		ret = min(ret, tmp);
	}
	return ret;
}

unsigned long query(unsigned int l, unsigned int r, unsigned int x){// x belong to s2, l&r belong to s
    // prume++;
    if (l > r) return 0;
    unsigned int mid = (l+r)/2;
	
	unsigned long resultTmp = 0;

    if (l == r) {
        //if (dist3(l, x) <= arcDis) {
        if (dist2(l, x) <= DIS) {
			return 1;
		}
        return 0;
    }
	//if (dist3(mid, x) <= arcDis){
	if (dist2(mid, x) <= DIS) {
		resultTmp += 1;
	}
	
    double distl = f(lc[mid], x), distr = f(rc[mid], x);
	
    if (distl <= arcDis){
        resultTmp += query(l, mid-1, x);
    }

    if (distr <= arcDis){
        resultTmp += query(mid+1, r, x);
    }
	return resultTmp;
}

void* threadQuery(void *arg){
	unsigned int l, r, x, xMax;
	l = ((ThreadArg*)arg)->l;
	r = ((ThreadArg*)arg)->r;
	x = ((ThreadArg*)arg)->x;
	xMax = ((ThreadArg*)arg)->xMax;

	unsigned long* result = (unsigned long*)malloc(sizeof(unsigned long));
	*result = 0;
	
	for (unsigned int i=x; i<=xMax; i++){
		*result += query(l, r, i);
	}
	
	pthread_exit(result);
}


struct LoadArg{
	unsigned int a, b;
	node* s;
	double* buf;
};

void* threadLoad(void* arg){
	unsigned int a, b;
	node* s;
	double* buf;    

	double ra, dec;

	s = ((LoadArg*)arg)->s;
	buf = ((LoadArg*)arg)->buf;
	a = ((LoadArg*)arg)->a;
	b = ((LoadArg*)arg)->b;

	for (unsigned int i = a; i<b; i++){
		ra = buf[i*2], dec = buf[i*2+1];
		//printf("%u %lf %lf\n", cnt, ra, dec);

		ra = ra*M_PI/180;
		dec = dec*M_PI/180;
        s[i+1].x = ra, s[i+1].y = dec;
	}

	pthread_exit(NULL);
}


void loadData(node* s, FILE* fp, unsigned int cntLimit, unsigned int & cnt){

	double* buf = (double*)malloc((cntLimit+1)*2*sizeof(double));
	cnt = fread(buf, sizeof(double), cntLimit*2, fp);
	cnt /= 2;
	//cnt = cntLimit;
	
	pthread_t thread[THREADNUM];
	LoadArg threadArg[THREADNUM];
	int res;
	void* threadRet;

	for (long i=0; i<curThreadNum; i++){
		threadArg[i].s = s;
		threadArg[i].buf = buf;
		threadArg[i].a = (unsigned int)(i*cnt/curThreadNum);
		threadArg[i].b = (unsigned int)((i+1)*cnt/curThreadNum);
	}


	for (unsigned int i=0; i<curThreadNum; i++){
		res = pthread_create(&thread[i], NULL, threadLoad, (void*)(&(threadArg[i])));
		if (res != 0){
			printf("Create thread failed");
			exit(res);
		}
	}


	for (unsigned int i=0; i<curThreadNum; i++){
		res = pthread_join(thread[i], &threadRet);
		if (res){
			printf("thread not joined.\n");
		}
	}


	free(buf);
	return;
}

int main(int argc, char* argv[]){

	time_t tStart = time(NULL), t;

    clock_t t1, t2, t3, t4;
    double ra, dec;

	int res;
	void* threadRet = NULL;


    char filePath[40] = "../data/short_wise"; // wise
    //char filePath[40] = "../data/short_psc_all";
    //char filePath[40] = "../data/short_sdssdr12";
    char filePath2[40] = "../data/short_sdssdr12";
    // char filePath2[40] = "../data/short_psc_all"; // 2MASS
    //char filePath2[40] = "../data/short_wise"; // wise
	// ra: 0-360, dec: -90-90

    unsigned int cnt = 0, cnt2 = 0, cntLimit = 747634026, cntLimit2 = 1231051050; //747634026
    //u 1231051050,  470992970, 747634026
    if (argc > 1) cntLimit = atoi(argv[1]);
    if (argc > 2) cntLimit2 = atoi(argv[2]);
    if (argc > 3) queryThreadNum = atoi(argv[3]);
	if (argc > 4) curThreadNum = atoi(argv[4]);
	if (argc > 5){
		if (strcmp(argv[5], "wise") == 0){
			strcpy(filePath, "../data/short_wise");
			cntLimit = 747634026;
		}else if(strcmp(argv[5], "sdss") == 0){
			strcpy(filePath, "../data/short_sdssdr12");
			cntLimit = 1231051050;
		}else{
			strcpy(filePath, "../data/short_psc_all");
			cntLimit = 470992970;
		}
	}

	if (argc > 6){
		if (strcmp(argv[6], "wise") == 0){
			strcpy(filePath2, "../data/short_wise");
			cntLimit2 = 747634026;
		}else if(strcmp(argv[6], "sdss") == 0){
			strcpy(filePath2, "../data/short_sdssdr12");
			cntLimit2 = 1231051050;
		}else{
			strcpy(filePath2, "../data/short_psc_all");
			cntLimit2 = 470992970;
		}
	}


    FILE* fp = fopen(filePath, "rb");
    FILE* fp2 = fopen(filePath2, "rb");
    if (fp == NULL){
        printf("[Error] No file.\n");
        return -1;
    }
	if (fp2 == NULL){
        printf("[Error] No file.\n");
        return -1;
    }

	// 分配空间
	t = time(NULL);
	s = (struct node*)malloc((cntLimit+1) * sizeof(node));
	s2 = (struct node*)malloc((cntLimit2+1) * sizeof(node));
	if (s2 == NULL || s == NULL){
		printf("[Error] Bad Malloc.\n");
	}

    //  load data
	printf("load\n");

	loadData(s, fp, cntLimit, cnt);
	loadData(s2, fp2, cntLimit2, cnt2);
	n = cnt;

	fclose(fp);
	fclose(fp2);
	
	// 之后再分配空间，节约资源
	cntLimit++;
	lc = (unsigned int*)malloc(cntLimit*sizeof(unsigned int));
	rc = (unsigned int*)malloc(cntLimit*sizeof(unsigned int));
	L = (double*)malloc(cntLimit*sizeof(double));
	R = (double*)malloc(cntLimit*sizeof(double));
	D = (double*)malloc(cntLimit*sizeof(double));
	U = (double*)malloc(cntLimit*sizeof(double));
	memset(lc, 0, cntLimit*sizeof(unsigned int));
	memset(rc, 0, cntLimit*sizeof(unsigned int));
	cntLimit--;

	double cost_t0 = time(NULL) - t;

	// build
	
	printf("build\n");
	t = time(NULL);
    
	BuildArg buildThreadMainArgs = {1, cnt, 0};
	pthread_t buildThreadMain;

	res = pthread_create(&buildThreadMain, NULL, threadBuild, (void*)(&(buildThreadMainArgs)));
	if (res != 0){
		printf("Create thread failed");
		exit(res);
	}
	res = pthread_join(buildThreadMain, &threadRet);
	if (res){
		printf("thread not joined.\n");
	}
	free(threadRet);
	
	double cost_t1 = time(NULL) - t;






	// query
	printf("query\n");
	t = time(NULL);
	pthread_t thread[THREADNUM];
	ThreadArg threadArg[THREADNUM];
	
	for (long i=0; i<queryThreadNum; i++){
		threadArg[i].l = 1;
		threadArg[i].r = cnt;
		threadArg[i].x = (unsigned int)(i*cnt2/queryThreadNum+1);
		threadArg[i].xMax = (unsigned int)((i+1)*cnt2/queryThreadNum);
	}


	
	
	for (unsigned int i=0; i<queryThreadNum; i++){
		// printf("creating %d thread.\n", i);
		res = pthread_create(&thread[i], NULL, threadQuery, (void*)(&(threadArg[i])));
		if (res != 0){
			printf("Create thread failed");
			exit(res);
		}
	}

	unsigned long queryResult = 0;

	for (unsigned int i=0; i<queryThreadNum; i++){
		res = pthread_join(thread[i], &threadRet);
		if (res){
			printf("thread not joined.\n");
		}
		queryResult += *((unsigned long*)(threadRet));
		free(threadRet);
	}
	double cost_t2 = time(NULL) - t;


	free(s);
	free(d);
	free(lc);
	free(rc);
	free(L);
	free(R);
	free(D);
	free(U);
	free(s2);





	double cost_t3 = time(NULL) - tStart;

    printf("[Info] Load time: %lf\n", cost_t0);
    printf("[Info] Build time: %lf\n", cost_t1);
    printf("[Info] Query time: %lf\n", cost_t2);
    printf("[Info] Total time: %lf\n", cost_t3);
    printf("[Info] Result: %ld\n", queryResult);
    cout << "[Info] Prume times: "<< prume << endl;

    return 0;
}



