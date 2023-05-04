#define _USE_MATH_DEFINES

#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define RAD M_PI/180

#define n1 40        	/*number of divisions in R direction*/
#define n2 40        	/*number of divisions in B direction*/
#define n3 40          	/*number of divisions in L direction*/

#define Lu 60.0         /*upper boundary for L*/
#define Ld 10.0         /*lower boundary for L*/
#define Bu 70.0 	/*upper boundary for B*/
#define Bd 20.0 	/*lower boundary for B*/

#define Ru 2.0 		/*upper boundary for R*/
#define Rd 1.0		/*lower boundary for R*/

#define GM 1.0
/*#define GM 3.986005e+14*/

/*parameters for SOR*/
#define omega 1.9
#define tol 1.0e-20     /*toleration*/
#define max_it 10000	 /*number of itterations*/

// na skusku subor:
// 1. Inicializacia prostredia
// 2. Jeden proces nacita subor
// 3. Poslu sa tie data
// 4. Iny proces ho vypise

int main()
{
	int nprocs = 0;
	int myRank = 0;

	// Initialize MPI environment
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs); // get num of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank); // assign ranks


	int i,j,k,it;
	double deltaR,deltaL,deltaB,deltaL1,deltaB1,res,pom,pom1,z,sigma,res3;
	double resLocal = 0, res3Local = 0;
	typedef double pole1[n1+2][n2+2][n3+2]; // + 2, aby sa dali zakomponovat aj okraje
	typedef double pole2[n1+2][n2+2][n3+2];
	static pole1 B,L,R;
	static pole2 T_B,T_L,T_R,s,u,deltag;
	static pole2 an,as,aw,ae,au,ad,ap,b,res2; // koeficienty v matici systemu "a"
	int meshSizeAll = (n1 + 2) * (n2 + 2) * (n3 + 2);
	int meshSizeSlice = (n2 + 2) * (n3 + 2); // = (n2 + 2) * (n3 + 2)

	FILE *fr, *fw;

	fw = fopen("stat.txt","w");

	deltaL=(Lu-Ld)/n3;
	deltaB=(Bu-Bd)/n2;
	deltaR=(Ru-Rd)/n1;
	deltaL1=deltaL*RAD*Rd;
	deltaB1=deltaB*RAD*Rd;
	//if (myRank == 0)
	//	printf("(dL1, dB1, dR) = (%.2lf, %.2lf, %.2lf)\n",deltaL1,deltaB1,deltaR);

	/*----------------------------------------------------------------------------*/
	/*vytvorenie 3D siete*/

	if (myRank == 0)
	{
		for (i = 1; i <= n1 + 1; i++) // + 1, aby sa dali spocitat stredy KO
			for (j = 1; j <= n2 + 1; j++)
				for (k = 1; k <= n3 + 1; k++)
				{
					L[i][j][k] = Ld + (k - 1) * deltaL;
					B[i][j][k] = Bd + (j - 1) * deltaB;
					R[i][j][k] = Rd + (i - 1) * deltaR;
				}
	}
	
	// Broadcast L, B, R to others
	MPI_Bcast(L, meshSizeAll, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(B, meshSizeAll, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(R, meshSizeAll, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/*----------------------------------------------------------------------------*/
	/*vypocet tazisk v BLR*/

	if (myRank == 0)
	{
		for (i = 1; i <= n1; i++)
			for (j = 1; j <= n2; j++)
				for (k = 1; k <= n3; k++)
				{
					T_L[i][j][k] = 0.5 * (L[i][j][k] + L[i][j][k + 1]);
					T_B[i][j][k] = 0.5 * (B[i][j][k] + B[i][j + 1][k]);
					T_R[i][j][k] = R[i][j][k] + 0.5 * deltaR;
				}

		/*----------------------------------------------------------------------------*/
		/*pridanie okrajovych tazisk*/

		for (i = 1; i <= n1; i++)
			for (j = 1; j <= n2; j++)
			{
				T_B[i][j][0] = T_B[i][j][1];
				T_B[i][j][n3 + 1] = T_B[i][j][n3];

				T_L[i][j][0] = T_L[i][j][1] - deltaL;
				T_L[i][j][n3 + 1] = T_L[i][j][n3] + deltaL;

				T_R[i][j][0] = T_R[i][j][1];
				T_R[i][j][n3 + 1] = T_R[i][j][n3];
			}

		for (i = 1; i <= n1; i++)
			for (k = 1; k <= n3; k++)
			{
				T_B[i][0][k] = T_B[i][1][k] - deltaB;
				T_B[i][n2 + 1][k] = T_B[i][n2][k] + deltaB;

				T_L[i][0][k] = T_L[i][1][k];
				T_L[i][n2 + 1][k] = T_L[i][n2][k];

				T_R[i][0][k] = T_R[i][1][k];
				T_R[i][n2 + 1][k] = T_R[i][n2][k];
			}

		// i -> radialny, j -> dozadu/dopredu, k -> vlavo/vpravo
		for (j = 1; j <= n2; j++)
			for (k = 1; k <= n3; k++)
			{
				T_B[0][j][k] = T_B[1][j][k]; // iba P0
				T_B[n1 + 1][j][k] = T_B[n1][j][k]; // iba Plast

				T_L[0][j][k] = T_L[1][j][k]; // iba P0
				T_L[n1 + 1][j][k] = T_L[n1][j][k]; // iba Plast

				T_R[0][j][k] = T_R[1][j][k] - deltaR;
				T_R[n1 + 1][j][k] = T_R[n1][j][k] + deltaR;
			}
	}

	// Broadcast T_L, T_B, T_R to others
	MPI_Bcast(T_B, meshSizeAll, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(T_L, meshSizeAll, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(T_R, meshSizeAll, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	

/*----------------------------------------------------------------------------*/
/*vynulovanie poli*/
// toto sa neparalelni
for (i=0;i<=n1+1;i++)
for (j=0;j<=n2+1;j++)
for (k=0;k<=n3+1;k++)
	{
	u[i][j][k] = 0.0;
	aw[i][j][k] = 0.0;
	ae[i][j][k] = 0.0;
	as[i][j][k] = 0.0;
	an[i][j][k] = 0.0;
	au[i][j][k] = 0.0;
	ad[i][j][k] = 0.0;
	ap[i][j][k] = 0.0;
	b[i][j][k] = 0.0;
	s[i][j][k] = 0.0;
	res2[i][j][k] = 0.0;
	}

/*----------------------------------------------------------------------------*/
/*vypocet indexov pre jednotilive procesy*/

int istart = 0, iend = 0;
int nlocal = 0, nlast = 0;
int inext = 0; // index of next process
int iprev = 0; // index of previous process

MPI_Status status;
MPI_Request req1, req2, req3, req4;

nlocal = (n1 + 2) / nprocs + 1;
nlast = (n1 + 2) - (nprocs - 1) * nlocal;

inext = myRank + 1;
iprev = myRank - 1;

if (myRank == 0)
	iprev = MPI_PROC_NULL;

if (myRank == nprocs - 1)
	inext = MPI_PROC_NULL;

istart = myRank * nlocal;

if (myRank == nprocs - 1)
	iend = istart + nlast - 1;
else
	iend = istart + nlocal - 1;


//// nlocal size variables
//nlocal = (n1 + 2) / nprocs + 1;
//nlast = (n1 + 2) - (nprocs - 1) * nlocal;
//
//// indices for next/previous process
//inext = myRank + 1;
//iprev = myRank - 1;
//
//if (myRank == 0) // P0 has no previous process
//	iprev = MPI_PROC_NULL;
//
//if (myRank == nprocs - 1) // P_last has no next process
//	inext = MPI_PROC_NULL;
//
//// istart/iend variables
//istart = myRank * nlocal + 1;
//
//// iend indices for all processes
//iend = istart + nlocal - 1;

// change istart for P0 and iend for P_last
//if (myRank == 0) // pre P0 chceme zacinat az od 1, nie od 0
	//istart = istart;

//if (myRank == nprocs - 1) // pre P_last chceme iend (n2 - 1) - 1
//{
//	iend = istart + nlast - 1;
//	//isendNext = iend + 1;
//	//irecvNext = iend + 1;
//}

// kontrolny vypis pre vsetky indexi
printf("P%d: istart: %d\tiend: %d\ninext: %d\tiprev: %d\n\n", 
	myRank, istart, iend, inext, iprev/*, isendNext, isendPrev, irecvNext, irecvPrev*/);

//MPI_Finalize();
//return 0;

/*----------------------------------------------------------------------------*/
/*vypocet koeficientov, tu netreba komunikovat medzi procesmi - kazdy si vypocita iba svoje*/
// pre P0 je istart = 1, pre P_last je iend = 20 (normalne by bolo 21, kedze velkost je 20 + 2)

if (myRank == 0)
{
	for (i = istart + 1; i <= iend; i++)
		for (j = 1; j <= n2; j++)
			for (k = 1; k <= n3; k++)
			{
				//Obsah W a E steny = deltaB/2*(Rspodne^2-Rvrchne^2), ... Obsah W a E steny sa meni len s R. RAD preraba stupne na radiany.
				// vzdialenosti medzi bodmi v smere W a E = deltaL*R*cos(B), ... vzdialenosti sa menia s R a B lebo cim sme blizsie k polom tak sa vzdianenosti zmensuju (lebo su kratsie rovnobezky)
				aw[i][j][k] = (deltaB / 2. * RAD * (R[i + 1][j][k] * R[i + 1][j][k] - R[i][j][k] * R[i][j][k])) / (deltaL * RAD * T_R[i][j][k] * cos(T_B[i][j][k] * RAD));
				ae[i][j][k] = aw[i][j][k];

				//Obsah N a S steny = deltaL/2*(Rspodne^2-Rvrchne^2)*cos(B), ... Obsah N a S steny sa meni s R ale aj s B lebo cim sme blizsie k polom tak sa obsahy stien zmensuju (lebo su kratsie rovnobezky)
				// vzdialenosti medzi bodmi v smere N a S = deltaB*R, ... vzdialenosti sa menia s R lebo poludniky su vzdy rovnako dlhe
				an[i][j][k] = (deltaL / 2. * RAD * (R[i + 1][j][k] * R[i + 1][j][k] - R[i][j][k] * R[i][j][k]) * cos(B[i][j + 1][k] * RAD)) / (deltaB * RAD * T_R[i][j][k]);
				as[i][j][k] = (deltaL / 2. * RAD * (R[i + 1][j][k] * R[i + 1][j][k] - R[i][j][k] * R[i][j][k]) * cos(B[i][j][k] * RAD)) / (deltaB * RAD * T_R[i][j][k]);

				//Obsah D steny = deltaL * (sin(Be) - sin(Bw)) * Rspodne^2
				//Obsah U steny = deltaL * (sin(Be) - sin(Bw)) * Rvrchne^2
				//vzdialenost medzi bodmi v smere U a D = deltaR
				ad[i][j][k] = deltaL * RAD * (sin(B[i + 1][j + 1][k] * RAD) - sin(B[i + 1][j][k] * RAD)) * R[i][j][k] * R[i][j][k] / deltaR;
				au[i][j][k] = deltaL * RAD * (sin(B[i][j + 1][k] * RAD) - sin(B[i][j][k] * RAD)) * R[i + 1][j][k] * R[i + 1][j][k] / deltaR;

				//vsetky aw, ae, an, as, ad, au by mali byt s minuskou ale v SOR sa to potom kompenzuje.
				ap[i][j][k] = aw[i][j][k] + ae[i][j][k] + as[i][j][k] + an[i][j][k] + au[i][j][k] + ad[i][j][k];
			}
}
else if (myRank == nprocs - 1)
{
	for (i = istart; i <= iend - 1; i++)
		for (j = 1; j <= n2; j++)
			for (k = 1; k <= n3; k++)
			{
				//Obsah W a E steny = deltaB/2*(Rspodne^2-Rvrchne^2), ... Obsah W a E steny sa meni len s R. RAD preraba stupne na radiany.
				// vzdialenosti medzi bodmi v smere W a E = deltaL*R*cos(B), ... vzdialenosti sa menia s R a B lebo cim sme blizsie k polom tak sa vzdianenosti zmensuju (lebo su kratsie rovnobezky)
				aw[i][j][k] = (deltaB / 2. * RAD * (R[i + 1][j][k] * R[i + 1][j][k] - R[i][j][k] * R[i][j][k])) / (deltaL * RAD * T_R[i][j][k] * cos(T_B[i][j][k] * RAD));
				ae[i][j][k] = aw[i][j][k];

				//Obsah N a S steny = deltaL/2*(Rspodne^2-Rvrchne^2)*cos(B), ... Obsah N a S steny sa meni s R ale aj s B lebo cim sme blizsie k polom tak sa obsahy stien zmensuju (lebo su kratsie rovnobezky)
				// vzdialenosti medzi bodmi v smere N a S = deltaB*R, ... vzdialenosti sa menia s R lebo poludniky su vzdy rovnako dlhe
				an[i][j][k] = (deltaL / 2. * RAD * (R[i + 1][j][k] * R[i + 1][j][k] - R[i][j][k] * R[i][j][k]) * cos(B[i][j + 1][k] * RAD)) / (deltaB * RAD * T_R[i][j][k]);
				as[i][j][k] = (deltaL / 2. * RAD * (R[i + 1][j][k] * R[i + 1][j][k] - R[i][j][k] * R[i][j][k]) * cos(B[i][j][k] * RAD)) / (deltaB * RAD * T_R[i][j][k]);

				//Obsah D steny = deltaL * (sin(Be) - sin(Bw)) * Rspodne^2
				//Obsah U steny = deltaL * (sin(Be) - sin(Bw)) * Rvrchne^2
				//vzdialenost medzi bodmi v smere U a D = deltaR
				ad[i][j][k] = deltaL * RAD * (sin(B[i + 1][j + 1][k] * RAD) - sin(B[i + 1][j][k] * RAD)) * R[i][j][k] * R[i][j][k] / deltaR;
				au[i][j][k] = deltaL * RAD * (sin(B[i][j + 1][k] * RAD) - sin(B[i][j][k] * RAD)) * R[i + 1][j][k] * R[i + 1][j][k] / deltaR;

				//vsetky aw, ae, an, as, ad, au by mali byt s minuskou ale v SOR sa to potom kompenzuje.
				ap[i][j][k] = aw[i][j][k] + ae[i][j][k] + as[i][j][k] + an[i][j][k] + au[i][j][k] + ad[i][j][k];
			}
}
else
{
	for (i = istart; i <= iend; i++)
		for (j = 1; j <= n2; j++)
			for (k = 1; k <= n3; k++)
			{
				//Obsah W a E steny = deltaB/2*(Rspodne^2-Rvrchne^2), ... Obsah W a E steny sa meni len s R. RAD preraba stupne na radiany.
				// vzdialenosti medzi bodmi v smere W a E = deltaL*R*cos(B), ... vzdialenosti sa menia s R a B lebo cim sme blizsie k polom tak sa vzdianenosti zmensuju (lebo su kratsie rovnobezky)
				aw[i][j][k] = (deltaB / 2. * RAD * (R[i + 1][j][k] * R[i + 1][j][k] - R[i][j][k] * R[i][j][k])) / (deltaL * RAD * T_R[i][j][k] * cos(T_B[i][j][k] * RAD));
				ae[i][j][k] = aw[i][j][k];

				//Obsah N a S steny = deltaL/2*(Rspodne^2-Rvrchne^2)*cos(B), ... Obsah N a S steny sa meni s R ale aj s B lebo cim sme blizsie k polom tak sa obsahy stien zmensuju (lebo su kratsie rovnobezky)
				// vzdialenosti medzi bodmi v smere N a S = deltaB*R, ... vzdialenosti sa menia s R lebo poludniky su vzdy rovnako dlhe
				an[i][j][k] = (deltaL / 2. * RAD * (R[i + 1][j][k] * R[i + 1][j][k] - R[i][j][k] * R[i][j][k]) * cos(B[i][j + 1][k] * RAD)) / (deltaB * RAD * T_R[i][j][k]);
				as[i][j][k] = (deltaL / 2. * RAD * (R[i + 1][j][k] * R[i + 1][j][k] - R[i][j][k] * R[i][j][k]) * cos(B[i][j][k] * RAD)) / (deltaB * RAD * T_R[i][j][k]);

				//Obsah D steny = deltaL * (sin(Be) - sin(Bw)) * Rspodne^2
				//Obsah U steny = deltaL * (sin(Be) - sin(Bw)) * Rvrchne^2
				//vzdialenost medzi bodmi v smere U a D = deltaR
				ad[i][j][k] = deltaL * RAD * (sin(B[i + 1][j + 1][k] * RAD) - sin(B[i + 1][j][k] * RAD)) * R[i][j][k] * R[i][j][k] / deltaR;
				au[i][j][k] = deltaL * RAD * (sin(B[i][j + 1][k] * RAD) - sin(B[i][j][k] * RAD)) * R[i + 1][j][k] * R[i + 1][j][k] / deltaR;

				//vsetky aw, ae, an, as, ad, au by mali byt s minuskou ale v SOR sa to potom kompenzuje.
				ap[i][j][k] = aw[i][j][k] + ae[i][j][k] + as[i][j][k] + an[i][j][k] + au[i][j][k] + ad[i][j][k];
			}
}


/*----------------------------------------------------------------------------*/
/*nacitanie okrajovych podmienok - Neumann*/

// presne riesenie je GM/R, derivacia v smere R je -GM/R^2
// na spodnej hranici je dana neumannova okrajova podmienka vsade inde dirichletova okrajova podmienka

if (myRank == 0)
{
	for (j = 1; j <= n2; j++) // iba proces P0
	{
		for (k = 1; k <= n3; k++)
		{
			//deltag je derivacia v radialnom smere (t.j. v smere normaly k D stene)
			deltag[1][j][k] = GM / (R[1][j][k] * R[1][j][k]);
			//kedze derivaciu v smere normaly k dolnej stene pozname, mozme ju prenasobit obsahom D steny a prehodit na pravu stranu
			b[1][j][k] = b[1][j][k] + deltaL * RAD * (sin(B[1][j + 1][k] * RAD) - sin(B[1][j][k] * RAD)) * R[1][j][k] * R[1][j][k] * deltag[1][j][k];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			ad[1][j][k] = 0;
			ap[1][j][k] = aw[1][j][k] + ae[1][j][k] + as[1][j][k] + an[1][j][k] + au[1][j][k];
		}
	}
}

/*----------------------------------------------------------------------------*/
/*nacitanie okrajovych podmienok - Dirichlet*/

if (myRank == nprocs - 1)
{
	for (j = 1; j <= n2; j++) // iba proces Plast
	{
		for (k = 1; k <= n3; k++)
		{
			//pom je presne riesenie v bode [n1+1][j][k]
			pom = GM / T_R[n1 + 1][j][k];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[n1][j][k] = b[n1][j][k] + pom * au[n1][j][k];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			au[n1][j][k] = 0;
		}
	}
}

if (myRank == 0)
{
	for (i = istart + 1; i <= iend; i++) // i = istart ... iend
		for (k = 1; k <= n3; k++)
		{
			//pom je presne riesenie v bode [i][0][k]
			pom = GM / T_R[i][0][k];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][1][k] = b[i][1][k] + pom * as[i][1][k];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			as[i][1][k] = 0;

			//pom je presne riesenie v bode [i][n2+1][k]
			pom = GM / T_R[i][n2 + 1][k];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][n2][k] = b[i][n2][k] + pom * an[i][n2][k];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			an[i][n2][k] = 0;
		}
}
else if (myRank == nprocs - 1)
{
	for (i = istart; i <= iend - 1; i++) // i = istart ... iend
		for (k = 1; k <= n3; k++)
		{
			//pom je presne riesenie v bode [i][0][k]
			pom = GM / T_R[i][0][k];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][1][k] = b[i][1][k] + pom * as[i][1][k];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			as[i][1][k] = 0;

			//pom je presne riesenie v bode [i][n2+1][k]
			pom = GM / T_R[i][n2 + 1][k];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][n2][k] = b[i][n2][k] + pom * an[i][n2][k];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			an[i][n2][k] = 0;
		}
}
else
{
	for (i = istart; i <= iend; i++) // i = istart ... iend
		for (k = 1; k <= n3; k++)
		{
			//pom je presne riesenie v bode [i][0][k]
			pom = GM / T_R[i][0][k];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][1][k] = b[i][1][k] + pom * as[i][1][k];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			as[i][1][k] = 0;

			//pom je presne riesenie v bode [i][n2+1][k]
			pom = GM / T_R[i][n2 + 1][k];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][n2][k] = b[i][n2][k] + pom * an[i][n2][k];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			an[i][n2][k] = 0;
		}
}

if (myRank == 0)
{
	for (i = istart + 1; i <= iend; i++) // i = istart ... iend
		for (j = 1; j <= n2; j++)
		{
			//pom je presne riesenie v bode [i][j][0]
			pom = GM / T_R[i][j][0];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][j][1] = b[i][j][1] + pom * aw[i][j][1];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			aw[i][j][1] = 0;

			//pom je presne riesenie v bode [i][j][n3+1]
			pom = GM / T_R[i][j][n3 + 1];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][j][n3] = b[i][j][n3] + pom * ae[i][j][n3];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			ae[i][j][n3] = 0;
		}
}
else if (myRank == nprocs - 1)
{
	for (i = istart; i <= iend - 1; i++) // i = istart ... iend
		for (j = 1; j <= n2; j++)
		{
			//pom je presne riesenie v bode [i][j][0]
			pom = GM / T_R[i][j][0];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][j][1] = b[i][j][1] + pom * aw[i][j][1];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			aw[i][j][1] = 0;

			//pom je presne riesenie v bode [i][j][n3+1]
			pom = GM / T_R[i][j][n3 + 1];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][j][n3] = b[i][j][n3] + pom * ae[i][j][n3];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			ae[i][j][n3] = 0;
		}
}
else
{
	for (i = istart; i <= iend; i++) // i = istart ... iend
		for (j = 1; j <= n2; j++)
		{
			//pom je presne riesenie v bode [i][j][0]
			pom = GM / T_R[i][j][0];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][j][1] = b[i][j][1] + pom * aw[i][j][1];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			aw[i][j][1] = 0;

			//pom je presne riesenie v bode [i][j][n3+1]
			pom = GM / T_R[i][j][n3 + 1];
			//kedze mame presne riesenie jedneho suseda, uz to nie je neznama a mozme ju prehodit na druhu stranu
			b[i][j][n3] = b[i][j][n3] + pom * ae[i][j][n3];
			//kedze sme ju prehodili na pravu stranu, treba ju nulovat medzi neznamymi
			ae[i][j][n3] = 0;
		}
}

/*----------------------------------------------------------------------------*/
/*riesenie sustavy rovnic*/
// RBSOR, tu uz sa bude musiet komunikovat medzi procesmi
pom=0.0;
it=0;

do
 {
	it=it+1;

	// UPDATE RED VOLUMES
	if (myRank == 0)
	{
		for (i = istart + 1; i <= iend; i++) // i = istart ... iend
			for (j = 1; j <= n2; j++)
				for (k = 1; k <= n3; k++)
				{
					if ((i + j + k) % 2 == 0) // kedze je tu Red-Black SOR, tak najprv sa ide cez Red konecne objemy
					{
						z = (b[i][j][k] + u[i][j][k + 1] * ae[i][j][k] +
							u[i][j][k - 1] * aw[i][j][k] +
							u[i][j + 1][k] * an[i][j][k] +
							u[i][j - 1][k] * as[i][j][k] +
							u[i + 1][j][k] * au[i][j][k] +
							u[i - 1][j][k] * ad[i][j][k]) / ap[i][j][k];

						u[i][j][k] = u[i][j][k] + omega * (z - u[i][j][k]);
					}
				}
	}
	else if (myRank == nprocs - 1)
	{
		for (i = istart; i <= iend - 1; i++) // i = istart ... iend
			for (j = 1; j <= n2; j++)
				for (k = 1; k <= n3; k++)
				{
					if ((i + j + k) % 2 == 0) // kedze je tu Red-Black SOR, tak najprv sa ide cez Red konecne objemy
					{
						z = (b[i][j][k] + u[i][j][k + 1] * ae[i][j][k] +
							u[i][j][k - 1] * aw[i][j][k] +
							u[i][j + 1][k] * an[i][j][k] +
							u[i][j - 1][k] * as[i][j][k] +
							u[i + 1][j][k] * au[i][j][k] +
							u[i - 1][j][k] * ad[i][j][k]) / ap[i][j][k];

						u[i][j][k] = u[i][j][k] + omega * (z - u[i][j][k]);
					}
				}
	}
	else
	{
		for (i = istart; i <= iend; i++) // i = istart ... iend
			for (j = 1; j <= n2; j++)
				for (k = 1; k <= n3; k++)
				{
					if ((i + j + k) % 2 == 0) // kedze je tu Red-Black SOR, tak najprv sa ide cez Red konecne objemy
					{
						z = (b[i][j][k] + u[i][j][k + 1] * ae[i][j][k] +
							u[i][j][k - 1] * aw[i][j][k] +
							u[i][j + 1][k] * an[i][j][k] +
							u[i][j - 1][k] * as[i][j][k] +
							u[i + 1][j][k] * au[i][j][k] +
							u[i - 1][j][k] * ad[i][j][k]) / ap[i][j][k];

						u[i][j][k] = u[i][j][k] + omega * (z - u[i][j][k]);
					}
				}
	}

	// komunikacia #1 po vypocte RED VOLUMES
	MPI_Isend(&u[iend], meshSizeSlice, MPI_DOUBLE, inext, 7, MPI_COMM_WORLD, &req1); // send up   (next P)
	MPI_Isend(&u[istart], meshSizeSlice, MPI_DOUBLE, iprev, 7, MPI_COMM_WORLD, &req2); // send down (prev P)
	MPI_Irecv(&u[iend + 1], meshSizeSlice, MPI_DOUBLE, inext, 7, MPI_COMM_WORLD, &req3); // recieve from up   (next P)
	MPI_Irecv(&u[istart - 1], meshSizeSlice, MPI_DOUBLE, iprev, 7, MPI_COMM_WORLD, &req4); // recieve from down (prev P)

	MPI_Wait(&req1, &status);
	MPI_Wait(&req2, &status);
	MPI_Wait(&req3, &status);
	MPI_Wait(&req4, &status);

	if (myRank == 0)
	{
		for (i = istart + 1; i <= iend; i++)
			for (j = 1; j <= n2; j++)
				for (k = 1; k <= n3; k++)
				{
					if ((i + j + k) % 2 == 1) // tu sa potom idu pocitat Black konecne objemy, pricom uz pozname updatenute vsetky Red konecne objemy
					{
						z = (b[i][j][k] + u[i][j][k + 1] * ae[i][j][k] +
							u[i][j][k - 1] * aw[i][j][k] +
							u[i][j + 1][k] * an[i][j][k] +
							u[i][j - 1][k] * as[i][j][k] +
							u[i + 1][j][k] * au[i][j][k] +
							u[i - 1][j][k] * ad[i][j][k]) / ap[i][j][k];

						u[i][j][k] = u[i][j][k] + omega * (z - u[i][j][k]);
					}
				}
	}
	else if (myRank == nprocs - 1)
	{
		for (i = istart; i <= iend - 1; i++)
			for (j = 1; j <= n2; j++)
				for (k = 1; k <= n3; k++)
				{
					if ((i + j + k) % 2 == 1) // tu sa potom idu pocitat Black konecne objemy, pricom uz pozname updatenute vsetky Red konecne objemy
					{
						z = (b[i][j][k] + u[i][j][k + 1] * ae[i][j][k] +
							u[i][j][k - 1] * aw[i][j][k] +
							u[i][j + 1][k] * an[i][j][k] +
							u[i][j - 1][k] * as[i][j][k] +
							u[i + 1][j][k] * au[i][j][k] +
							u[i - 1][j][k] * ad[i][j][k]) / ap[i][j][k];

						u[i][j][k] = u[i][j][k] + omega * (z - u[i][j][k]);
					}
				}
	}
	else
	{
		for (i = istart; i <= iend; i++)
			for (j = 1; j <= n2; j++)
				for (k = 1; k <= n3; k++)
				{
					if ((i + j + k) % 2 == 1) // tu sa potom idu pocitat Black konecne objemy, pricom uz pozname updatenute vsetky Red konecne objemy
					{
						z = (b[i][j][k] + u[i][j][k + 1] * ae[i][j][k] +
							u[i][j][k - 1] * aw[i][j][k] +
							u[i][j + 1][k] * an[i][j][k] +
							u[i][j - 1][k] * as[i][j][k] +
							u[i + 1][j][k] * au[i][j][k] +
							u[i - 1][j][k] * ad[i][j][k]) / ap[i][j][k];

						u[i][j][k] = u[i][j][k] + omega * (z - u[i][j][k]);
					}
				}
	}
		

	// komunikacia #2 po vypocte BLACK VOLUMES
	MPI_Isend(&u[iend], meshSizeSlice, MPI_DOUBLE, inext, 7, MPI_COMM_WORLD, &req1); // send up   (next P)
	MPI_Isend(&u[istart], meshSizeSlice, MPI_DOUBLE, iprev, 7, MPI_COMM_WORLD, &req2); // send down (prev P)
	MPI_Irecv(&u[iend + 1], meshSizeSlice, MPI_DOUBLE, inext, 7, MPI_COMM_WORLD, &req3); // recieve from up   (next P)
	MPI_Irecv(&u[istart - 1], meshSizeSlice, MPI_DOUBLE, iprev, 7, MPI_COMM_WORLD, &req4); // recieve from down (prev P)

	MPI_Wait(&req1, &status);
	MPI_Wait(&req2, &status);
	MPI_Wait(&req3, &status);
	MPI_Wait(&req4, &status);

	res = res3 = 0.0;
	resLocal = 0.0;
	res3Local = 0.0;

	// compute residuals
	if (myRank == 0)
	{
		for (i = istart + 1; i <= iend; i++)
			for (j = 1; j <= n2; j++)
				for (k = 1; k <= n3; k++)
				{
					pom = (u[i][j][k] * ap[i][j][k] - u[i][j][k + 1] * ae[i][j][k] -
						u[i][j][k - 1] * aw[i][j][k] -
						u[i][j + 1][k] * an[i][j][k] -
						u[i][j - 1][k] * as[i][j][k] -
						u[i + 1][j][k] * au[i][j][k] -
						u[i - 1][j][k] * ad[i][j][k] - b[i][j][k]);

					resLocal = resLocal + pom * pom;

					pom1 = GM / T_R[i][j][k] - u[i][j][k];
					res3Local = res3Local + pom1 * pom1;
				}
	}
	else if (myRank == nprocs - 1)
	{
		for (i = istart; i <= iend - 1; i++)
			for (j = 1; j <= n2; j++)
				for (k = 1; k <= n3; k++)
				{
					pom = (u[i][j][k] * ap[i][j][k] - u[i][j][k + 1] * ae[i][j][k] -
						u[i][j][k - 1] * aw[i][j][k] -
						u[i][j + 1][k] * an[i][j][k] -
						u[i][j - 1][k] * as[i][j][k] -
						u[i + 1][j][k] * au[i][j][k] -
						u[i - 1][j][k] * ad[i][j][k] - b[i][j][k]);

					resLocal = resLocal + pom * pom;

					pom1 = GM / T_R[i][j][k] - u[i][j][k];
					res3Local = res3Local + pom1 * pom1;
				}
	}
	else
	{
		for (i = istart; i <= iend; i++)
			for (j = 1; j <= n2; j++)
				for (k = 1; k <= n3; k++)
				{
					pom = (u[i][j][k] * ap[i][j][k] - u[i][j][k + 1] * ae[i][j][k] -
						u[i][j][k - 1] * aw[i][j][k] -
						u[i][j + 1][k] * an[i][j][k] -
						u[i][j - 1][k] * as[i][j][k] -
						u[i + 1][j][k] * au[i][j][k] -
						u[i - 1][j][k] * ad[i][j][k] - b[i][j][k]);

					resLocal = resLocal + pom * pom;

					pom1 = GM / T_R[i][j][k] - u[i][j][k];
					res3Local = res3Local + pom1 * pom1;
				}
	}

	// Reduce for resiudals variables res and res3
	MPI_Allreduce(&resLocal,  &res,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&res3Local, &res3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	if (myRank == 0) // only master prints message
	{
		res3 = sqrt(res3 / (n1 * n2 * n3));
		printf("\t\t%d\t%.12lf\n", it, res);
	}

} while ((res > tol) && (it < max_it));


/*----------------------------------------------------------------------------*/
/*porovnanie s presnym riesenim a zapis vysledkov*/

sigma=0.0;
pom=0.0;

double pomLocal = 0.0;

if (myRank == 0)
{
	for (i = istart + 1; i <= iend; i++)
		for (j = 1; j <= n2; j++)
			for (k = 1; k <= n3; k++)
			{
				res2[i][j][k] = (GM / T_R[i][j][k]) - u[i][j][k];
				pomLocal = pomLocal + res2[i][j][k] * res2[i][j][k];
				//fprintf(fw,"%d %d %d\t%.7lf\t%.9lf\n",i,j,k,u[i][j][k],res2[i][j][k]);
			}
}
else if (myRank == nprocs - 1)
{
	for (i = istart; i <= iend - 1; i++)
		for (j = 1; j <= n2; j++)
			for (k = 1; k <= n3; k++)
			{
				res2[i][j][k] = (GM / T_R[i][j][k]) - u[i][j][k];
				pomLocal = pomLocal + res2[i][j][k] * res2[i][j][k];
				//fprintf(fw,"%d %d %d\t%.7lf\t%.9lf\n",i,j,k,u[i][j][k],res2[i][j][k]);
			}
}
else
{
	for (i = istart; i <= iend; i++)
		for (j = 1; j <= n2; j++)
			for (k = 1; k <= n3; k++)
			{
				res2[i][j][k] = (GM / T_R[i][j][k]) - u[i][j][k];
				pomLocal = pomLocal + res2[i][j][k] * res2[i][j][k];
				//fprintf(fw,"%d %d %d\t%.7lf\t%.9lf\n",i,j,k,u[i][j][k],res2[i][j][k]);
			}
}

MPI_Allreduce(&pomLocal, &pom, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

if (myRank == 0)
{
	sigma = sqrt(pom / (n1 * n2 * n3)); // kazdy svoju lokalnu sigmu, potom MPI_Reduce cez MPI_SUM do globalnej sigmy
	printf("sigma = %.20lf\n", sigma);

}


// len P0, lebo sa vypisuje prva vrstva
if (myRank == 0)
{
	i = 1;
	for (j = 1; j <= n2; j++)
		for (k = 1; k <= n3; k++)
		{
			res2[i][j][1] = (GM / T_R[i][j][k]) - u[i][j][k];
			fprintf(fw, "%d %d %d\t%.7lf\t%.9lf\t%.7lf\n", i, j, k, u[i][j][k], res2[i][j][k], (GM / T_R[i][j][k]));
		}
}
/*
i=n1/2;
     for (j=1;j<=n2;j++)
       for (k=1;k<=n3;k++)
        {
	res2[i][j][1]=(GM/T_R[i][j][k])-u[i][j][k];
	fprintf(fw,"%d %d %d\t%.7lf\t%.9lf\t%.7lf\n",i,j,k,u[i][j][j],res2[i][j][j],(GM/T_R[i][j][k]));
	}

i=n1;
     for (j=1;j<=n2;j++)
       for (k=1;k<=n3;k++)
        {
	res2[i][j][1]=(GM/T_R[i][j][k])-u[i][j][k];
	fprintf(fw,"%d %d %d\t%.7lf\t%.9lf\t%.7lf\n",i,j,k,u[i][j][j],res2[i][j][j],(GM/T_R[i][j][k]));
	}
*/
fclose(fw);

MPI_Finalize();

return 0;
}
