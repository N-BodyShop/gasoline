#ifndef SSDEFS_HINCLUDED
#define SSDEFS_HINCLUDED

#ifndef N_DIM
#define N_DIM 3
#endif

struct ss_head {
	double	time;
	int		n_data;
	int		pad;	/* unused: pad to 8 byte boundary */
	};

struct ss_data {
	double	mass;
	double	radius;
	double	pos[N_DIM];
	double	vel[N_DIM];
	double	spin[N_DIM];
	int		color;
	int		pad;	/* unused: pad to 8 byte boundary */
	};

/* Object color identifiers */

#define SUN				5	/* Yellow */
#define JUPITER			2	/* Red */
#define SATURN			11	/* Khaki */
#define URANUS			6	/* Magenta */
#define NEPTUNE			7	/* Cyan */
#define PLANETESIMAL	3	/* Green */

#endif /* !SSDEFS_HINCLUDED */
