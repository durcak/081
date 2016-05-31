#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define BODIES 5
/* Jupiter, Io, Europa, Ganymed, Callisto */

/* hmotnost */
const double mass[] = {
	1.9e27, 8.9e22, 4.8e22, 1.48e23, 1.08e23
};

/* polomer planety/mesice */
const double radius[] = {
	72e6, 1.82e6, 1.57e6, 2.63e6, 2.41e6
};

/* polomer obezne drahy */
const double axis[] = {
	422e6, 671e6, 1070e6, 1882e6
};

/* doba obehu */
const double period[] = {
	152928.00, 306720.00, 617760.00, 1442880.0
};

/* gravitacni konstanta */
const double kappa = 6.67e-11;

/* pocatecni poloha ulomku (x_1), x_2 = 0 */
#define START 2e9

/* vzdalenost, kterou povazujeme za unik */
#define LEAVE 10e9

/* pocatecni rychlost */
#define START_V 7000.
#define START_ALPHA (M_PI*.7)

#define SEED 123456

/* pocet ulomku */
#define COUNT_MAX	1000
int	count = 1;

/* casovy limit */
double limit = 1e6;

/* velikost kroku integrace */
double delta_t = 1;

init_state(double x[][2], double v[][2])
{
	double	vx = cos(START_ALPHA)*START_V,
		vy = sin(START_ALPHA)*START_V;

	int	i;

	fprintf(stderr,"init\n");
	for (i=0; i<count; i++) {
		x[i][0] = START;
		x[i][1] = 0;

	/* nahodny vektor, smer uniformni, velikost exponencialne klesajici, max 1/10
	   pocatecni rychlosti komety */
		double vo = 2.*M_PI*drand48(),
		      vr = START_V*1.0*exp(-15.*drand48());

		v[i][0] = vx + vr*cos(vo);
		v[i][1] = vy + vr*sin(vo);

		fprintf(stderr,"%5d: %f %f %f %f\n",i,x[i][0],x[i][1],v[i][0],v[i][1]);
	}
	fprintf(stderr,"end init\n");
}

/* vyhodnoceni pohybove rovnice, tj. vypocet casove derivace */
eval_f(
	double x[2],	/* aktualni poloha ulomku */
	double v[2],	/* aktualni rychlost */
	double m[][2],	/* pozice mesicu */
	double dx[2],		/* vypoctena derivace polohy */
	double dv[2])		/* vypoctena derivace rychlosti */
{
	dx[0] = v[0];
	dx[1] = v[1];

	/* kvadrat a reciproka vzdalenost od Jupitera */
	double r2 = x[0]*x[0] + x[1]*x[1],rr = 1./sqrt(r2);	

	double a = kappa * mass[0]/r2;	/* gravitacni zrychleni v danem miste */

	/* odpovidajici vektor */
	dv[0] = -x[0]*rr * a;
	dv[1] = -x[1]*rr * a;

	int	j;
	for (j=0; j<BODIES-1; j++) {
		double	d[2],d2,rd;

		/* vektor od ulomku k mesici j */
		d[0] = m[j][0] - x[0];
		d[1] = m[j][1] - x[1];

		/* kvadrat a reciproka velikost */
		d2 = d[0]*d[0] + d[1]*d[1];
		rd = 1./sqrt(d2);

		double a = kappa * mass[j+1]/d2;	/* gravitacni zrychleni mesice j */

		d[0] *= rd * a;	/* odpovidajici vektor */
		d[1] *= rd * a;

		dv[0] += d[0];	/* kumulace pres vsechny mesice */
		dv[1] += d[1];
	}
}

int check_active(int i, double t, double x[2], double v[2], double m[][2])
{
	double	r2 = x[0]*x[0] + x[1]*x[1];
	int	j;

	/* vyletel z dosahu */
	if (r2 > LEAVE*LEAVE) {
		double	d = sqrt(r2) - LEAVE, /* jak jsme daleko */
			vv = sqrt(v[0]*v[0] + v[1]*v[1]); /* jakou rychlosti */
		
		printf ("out %d: -1 %18.3f %13.7g %13.7g\n",i,
				t - d/vv,
				x[0] - v[0]*d/vv, x[1] - v[1]*d/vv);
		return 0;
	}
	else for (j=0; j<BODIES; j++) {
		double	d[2],d2;

		if (j == 0) {	/* vektor od ulomku k Jupiteru */
			d[0] = x[0];
			d[1] = x[1];
		}
		else { /* vektor od ulomku k mesici j */
			d[0] = x[0] - m[j-1][0];
			d[1] = x[1] - m[j-1][1];
		}

		d2 = d[0]*d[0] + d[1]*d[1];
		if (d2 < radius[j]*radius[j]) { /* trefil se */
		/* parametr p -- jak daleko se vydat z x po v, abychom se vratili na kruznici.
		 * vyjde z toho kvadraticka rovnice */
			double	a = v[0]*v[0] + v[1]*v[1],
				b = 2.*(v[0]*d[0] + v[1]*d[1]),
				c = d2 - radius[j]*radius[j],
				det = b*b - 4.*a*c,
				p = (-b - sqrt(det))/(2.*a);
			
			printf("out %d: %d %18.3f %13.7g %13.7g\n",i,j,
					t + p,
					x[0] + v[0]*p,
					x[1] + v[1]*p);
			return 0;
		}
	}
	return 1;
}

void usage(int argc, char **argv)
{
	fprintf(stderr,"usage: %s [options]\n"
			" -s time-step\n"
			" -t state-trace-modulo\n"
			" -n number-of-pieces\n"
			" -l time-limit\n"
			" -m method (0 euler, 1 leapfrog, 2 runge-kutta)\n",argv[0]);
}

int	trace = 0;

main(int argc, char **argv)
{
	double	t = 0;	/* globalni cas */

	double	x[COUNT_MAX][2], v[COUNT_MAX][2];

	int	nactive, active[COUNT_MAX];		/* ktere a kolik ulomku se jeste hybe */
	int	i;
	long	step = 0;

	int	opt,method = 0;

	printf("options: ");
	for (i=0; i<argc; i++) printf("%s ",argv[i]);
	printf("\n");

	while ((opt = getopt(argc,argv,"l:s:t:n:m:")) != EOF) switch (opt) {
		case 's': delta_t = atof(optarg); break;
		case 't': trace = atoi(optarg); break;
		case 'n': count = atoi(optarg); break;
		case 'm': method = atoi(optarg); break;
		case 'l': limit = atof(optarg); break;
		default: usage(argc,argv); exit(1);
	}

	srand48(SEED);
	nactive = count;
	for (i=0; i<count; i++) active[i] = 1;
	init_state(x,v);

	while (nactive > 0 && t<limit) {
		double	omega[BODIES-1],m[BODIES-1][2];

		for (i=0; i<BODIES-1; i++) {
			omega[i] = 2.*M_PI*t/period[i];
			m[i][0] = axis[i]*cos(omega[i]);
			m[i][1] = axis[i]*sin(omega[i]);
		}

		for (i=0; i<count; i++) if (active[i]) {
			double	dx[2],dv[2],dv2[2];
			double	k[5][4],y[5][4];
			int	j;

			switch (method) {
				case 0:
		/* krok dopredne Eulerovy integrace */
					eval_f(x[i],v[i],m,dx,dv);
					x[i][0] += dx[0]*delta_t;
					x[i][1] += dx[1]*delta_t;
					v[i][0] += dv[0]*delta_t;
					v[i][1] += dv[1]*delta_t;
					break;
				case 1: /*leapfrog */
					break;
				case 2: /* runge-kutta */

     					eval_f(x[i],v[i],m,dx,dv);
					k[0][1] = delta_t * dx[0];
					k[1][1] = delta_t * dx[1];

					pom[0]  = x[i][0] + delta_t/2;
					pom[1]  = x[i][1] + delta_t/2;
					pom2[0] = v[i][0] + k[0][1]/2;
					pom2[1] = v[i][1] + k[1][1]/2; 
				  	eval_f(pom, pom2, m, dx, dv);

					k[0][2] = delta_t * dx[0];
					k[1][2] = delta_t * dx[1];

					pom[0]  = x[i][0] + delta_t/2;
					pom[1]  = x[i][1] + delta_t/2;
					pom2[0] = v[i][0] + k[0][2]/2;
					pom2[1] = v[i][1] + k[1][2]/2; 
				  	eval_f(pom, pom2, m, dx, dv);

					k[0][3] = delta_t * dx[0];
					k[1][3] = delta_t * dx[1];

					pom[0]  = x[i][0] + delta_t;
					pom[1]  = x[i][1] + delta_t;
					pom2[0] = v[i][0] + k[0][2];
					pom2[1] = v[i][1] + k[1][2]; 
				  	eval_f(pom, pom2, m, dx, dv);

					k[0][4] = delta_t * dx[0];
					k[1][4] = delta_t * dx[1];
		
					dx[0] = (k[0][1] + 2.0 * k[0][2] + 2.0 * k[0][3] + k[0][4])/6.0;
					dx[1] = (k[1][1] + 2.0 * k[1][2] + 2.0 * k[1][3] + k[1][4])/6.0;
					x[i][0] += dx[0];
					x[i][1] += dx[1];

					//dv[0] = (k[2][1] + 2.0 * k[2][2] + 2.0 * k[2][3] + k[2][4])/6.0;
					//dv[1] = (k[3][1] + 2.0 * k[3][2] + 2.0 * k[3][3] + k[3][4])/6.0;
					v[i][0] += dv[0];
					v[i][1] += dv[1];

				default:
					fprintf(stderr,"unknown method\n"); exit(1);
			}

			if (!check_active(i,t,x[i],v[i],m)) {
				active[i] = 0;
				nactive--;
			}
		}

		if (trace && ((step % trace) == 0)) {
			printf ("stat: ");
			for (i=0; i<BODIES-1; i++) printf("%f %f ",m[i][0],m[i][1]);
			for (i=0; i<count; i++)
				if (active[i]) printf("%f %f ",x[i][0],x[i][1]);
				else printf("0 0 ");
			printf("\n");
		}
	
		t += delta_t;
		step++;
	}
	for (i=0; i<count; i++)
		printf("final: %c %5d: %f %f %f %f\n",active[i] ? '*' : ' ',i,x[i][0],x[i][1],v[i][0],v[i][1]);
}
