#ifndef SSDEFS_HINCLUDED
#define SSDEFS_HINCLUDED

#ifndef N_DIM
#define N_DIM 3
#endif

/*
 ** Following structures intended for use with xdr I/O routines. Note that
 ** internal representations may differ from format written to disk. As a
 ** result, for now use hardwired sizes for seeking. Is there an xdr routine
 ** to figure this out automatically?...
 */

struct ss_head {
	double	time;
	int		n_data;
	int		pad;	/* unused: pad to 8 byte boundary */
	};

#define SS_HEAD_SIZE 16	/* sizeof(ss_head) */

struct ss_data {
	double	mass;
	double	radius;
	double	pos[N_DIM];
	double	vel[N_DIM];
	double	spin[N_DIM];
	int		color;
	int		pad;	/* unused: pad to 8 byte boundary */
	};

#define SS_DATA_SIZE 96	/* sizeof(ss_data) */

/* Object color identifiers */

#define SUN				5	/* Yellow */
#define JUPITER			2	/* Red */
#define SATURN			11	/* Khaki */
#define URANUS			6	/* Magenta */
#define NEPTUNE			7	/* Cyan */
#define PLANETESIMAL	3	/* Green */

#endif /* !SSDEFS_HINCLUDED */
