#if 0
	shc Version 3.9.6, Generic Shell Script Compiler
	GNU GPL Version 3 Md Jahidul Hamid <jahidulhamid@yahoo.com>

	shc -f bamplot_sh_wrap.sh 
#endif

static  char data [] = 
#define      xecc_z	15
#define      xecc	((&data[0]))
	"\210\002\237\337\041\273\337\157\202\305\121\060\156\130\114\040"
	"\143"
#define      chk1_z	22
#define      chk1	((&data[18]))
	"\313\044\223\113\311\221\117\243\256\077\241\307\037\266\101\017"
	"\075\031\173\223\163\313\241\142\313\276\262"
#define      lsto_z	1
#define      lsto	((&data[44]))
	"\006"
#define      opts_z	1
#define      opts	((&data[45]))
	"\155"
#define      text_z	25
#define      text	((&data[50]))
	"\207\004\074\220\032\031\054\334\001\030\243\226\315\100\221\114"
	"\273\371\150\115\310\107\037\076\243\233\075\253\173\107\270"
#define      msg2_z	19
#define      msg2	((&data[80]))
	"\363\344\277\116\341\046\164\343\263\333\157\011\337\323\063\054"
	"\332\166\003\113\002\113"
#define      msg1_z	65
#define      msg1	((&data[103]))
	"\143\034\041\305\217\042\314\130\151\324\162\260\340\060\277\021"
	"\174\266\141\310\130\377\322\006\326\315\212\157\023\207\233\103"
	"\140\131\073\146\216\302\245\353\316\254\155\345\056\044\150\173"
	"\164\231\067\360\331\205\163\015\374\340\054\154\037\076\031\175"
	"\242\044\120\232\204\347\337\170\100\207"
#define      tst1_z	22
#define      tst1	((&data[173]))
	"\112\273\071\022\175\075\351\076\313\074\275\161\223\322\212\335"
	"\155\277\050\114\032\374"
#define      inlo_z	3
#define      inlo	((&data[195]))
	"\240\043\260"
#define      rlax_z	1
#define      rlax	((&data[198]))
	"\244"
#define      chk2_z	19
#define      chk2	((&data[203]))
	"\257\207\317\314\167\055\132\160\027\122\162\314\275\376\157\111"
	"\333\376\233\122\326\372\065"
#define      tst2_z	19
#define      tst2	((&data[225]))
	"\216\374\341\274\312\023\067\262\307\200\331\001\155\225\316\071"
	"\047\136\201\262\014\117\362\030"
#define      pswd_z	256
#define      pswd	((&data[292]))
	"\342\060\100\152\270\112\205\007\051\101\302\220\001\272\207\010"
	"\152\017\330\066\212\214\302\031\210\244\013\240\247\303\241\211"
	"\363\342\364\254\054\172\264\125\273\167\345\274\062\155\325\334"
	"\025\250\107\266\064\351\304\034\337\100\212\051\026\124\044\365"
	"\251\013\064\131\005\060\076\103\362\106\250\273\116\176\230\143"
	"\047\337\032\133\310\336\170\250\037\002\322\065\127\366\053\000"
	"\002\140\132\010\220\230\113\202\337\364\076\055\163\326\221\232"
	"\266\253\366\177\212\156\047\252\161\371\337\310\360\013\310\363"
	"\153\042\373\373\273\106\176\232\073\275\310\256\223\131\110\111"
	"\005\076\311\217\254\360\071\035\352\031\346\333\045\256\316\220"
	"\321\311\214\215\020\013\050\113\310\360\371\133\111\101\245\117"
	"\200\156\336\054\137\030\112\112\062\061\045\127\340\364\347\261"
	"\275\164\077\315\177\147\031\107\127\022\243\241\124\111\360\324"
	"\270\317\001\030\350\114\142\032\175\210\161\135\175\131\017\072"
	"\316\116\010\115\265\041\225\015\064\071\256\210\203\237\135\074"
	"\157\137\124\127\253\267\162\050\100\344\206\275\075\225\367\014"
	"\344\000\131\232\042\357\247\126\051\126\337\254\365\075\350\145"
	"\234\075\274\107\364\056\160\065\022\366\362\120\214\352\134\161"
	"\352\266\013\015\246\262\143\317\010\103\174\376\200\144"
#define      date_z	1
#define      date	((&data[548]))
	"\332"
#define      shll_z	10
#define      shll	((&data[551]))
	"\052\225\042\103\210\052\377\353\073\366\253\276"/* End of data[] */;
#define      hide_z	4096
#define DEBUGEXEC	0	/* Define as 1 to debug execvp calls */
#define TRACEABLE	1	/* Define as 1 to enable ptrace the executable */
#define BUSYBOXON	0	/* Define as 1 to enable work with busybox */

/* rtc.c */

#include <sys/stat.h>
#include <sys/types.h>

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* 'Alleged RC4' */

static unsigned char stte[256], indx, jndx, kndx;

/*
 * Reset arc4 stte. 
 */
void stte_0(void)
{
	indx = jndx = kndx = 0;
	do {
		stte[indx] = indx;
	} while (++indx);
}

/*
 * Set key. Can be used more than once. 
 */
void key(void * str, int len)
{
	unsigned char tmp, * ptr = (unsigned char *)str;
	while (len > 0) {
		do {
			tmp = stte[indx];
			kndx += tmp;
			kndx += ptr[(int)indx % len];
			stte[indx] = stte[kndx];
			stte[kndx] = tmp;
		} while (++indx);
		ptr += 256;
		len -= 256;
	}
}

/*
 * Crypt data. 
 */
void arc4(void * str, int len)
{
	unsigned char tmp, * ptr = (unsigned char *)str;
	while (len > 0) {
		indx++;
		tmp = stte[indx];
		jndx += tmp;
		stte[indx] = stte[jndx];
		stte[jndx] = tmp;
		tmp += stte[indx];
		*ptr ^= stte[tmp];
		ptr++;
		len--;
	}
}

/* End of ARC4 */

/*
 * Key with file invariants. 
 */
int key_with_file(char * file)
{
	struct stat statf[1];
	struct stat control[1];

	if (stat(file, statf) < 0)
		return -1;

	/* Turn on stable fields */
	memset(control, 0, sizeof(control));
	control->st_ino = statf->st_ino;
	control->st_dev = statf->st_dev;
	control->st_rdev = statf->st_rdev;
	control->st_uid = statf->st_uid;
	control->st_gid = statf->st_gid;
	control->st_size = statf->st_size;
	control->st_mtime = statf->st_mtime;
	control->st_ctime = statf->st_ctime;
	key(control, sizeof(control));
	return 0;
}

#if DEBUGEXEC
void debugexec(char * sh11, int argc, char ** argv)
{
	int i;
	fprintf(stderr, "shll=%s\n", sh11 ? sh11 : "<null>");
	fprintf(stderr, "argc=%d\n", argc);
	if (!argv) {
		fprintf(stderr, "argv=<null>\n");
	} else { 
		for (i = 0; i <= argc ; i++)
			fprintf(stderr, "argv[%d]=%.60s\n", i, argv[i] ? argv[i] : "<null>");
	}
}
#endif /* DEBUGEXEC */

void rmarg(char ** argv, char * arg)
{
	for (; argv && *argv && *argv != arg; argv++);
	for (; argv && *argv; argv++)
		*argv = argv[1];
}

void chkenv_end(void);

int chkenv(int argc)
{
	char buff[512];
	unsigned long mask, m;
	int l, a, c;
	char * string;
	extern char ** environ;

	mask = (unsigned long)getpid();
	stte_0();
	 key(&chkenv, (void*)&chkenv_end - (void*)&chkenv);
	 key(&data, sizeof(data));
	 key(&mask, sizeof(mask));
	arc4(&mask, sizeof(mask));
	sprintf(buff, "x%lx", mask);
	string = getenv(buff);
#if DEBUGEXEC
	fprintf(stderr, "getenv(%s)=%s\n", buff, string ? string : "<null>");
#endif
	l = strlen(buff);
	if (!string) {
		/* 1st */
		sprintf(&buff[l], "=%lu %d", mask, argc);
		putenv(strdup(buff));
		return 0;
	}
	c = sscanf(string, "%lu %d%c", &m, &a, buff);
	if (c == 2 && m == mask) {
		/* 3rd */
		rmarg(environ, &string[-l - 1]);
		return 1 + (argc - a);
	}
	return -1;
}

void chkenv_end(void){}

#if !TRACEABLE

#define _LINUX_SOURCE_COMPAT
#include <sys/ptrace.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <signal.h>
#include <stdio.h>
#include <unistd.h>

#if !defined(PTRACE_ATTACH) && defined(PT_ATTACH)
#	define PTRACE_ATTACH	PT_ATTACH
#endif
void untraceable(char * argv0)
{
	char proc[80];
	int pid, mine;

	switch(pid = fork()) {
	case  0:
		pid = getppid();
		/* For problematic SunOS ptrace */
#if defined(__FreeBSD__)
		sprintf(proc, "/proc/%d/mem", (int)pid);
#else
		sprintf(proc, "/proc/%d/as",  (int)pid);
#endif
		close(0);
		mine = !open(proc, O_RDWR|O_EXCL);
		if (!mine && errno != EBUSY)
			mine = !ptrace(PTRACE_ATTACH, pid, 0, 0);
		if (mine) {
			kill(pid, SIGCONT);
		} else {
			perror(argv0);
			kill(pid, SIGKILL);
		}
		_exit(mine);
	case -1:
		break;
	default:
		if (pid == waitpid(pid, 0, 0))
			return;
	}
	perror(argv0);
	_exit(1);
}
#endif /* !TRACEABLE */

char * xsh(int argc, char ** argv)
{
	char * scrpt;
	int ret, i, j;
	char ** varg;
	char * me = argv[0];
	if (me == NULL) { me = getenv("_"); }
	if (me == 0) { fprintf(stderr, "E: neither argv[0] nor $_ works."); exit(1); }

	ret = chkenv(argc);
	stte_0();
	 key(pswd, pswd_z);
	arc4(msg1, msg1_z);
	arc4(date, date_z);
	if (date[0] && (atoll(date)<time(NULL)))
		return msg1;
	arc4(shll, shll_z);
	arc4(inlo, inlo_z);
	arc4(xecc, xecc_z);
	arc4(lsto, lsto_z);
	arc4(tst1, tst1_z);
	 key(tst1, tst1_z);
	arc4(chk1, chk1_z);
	if ((chk1_z != tst1_z) || memcmp(tst1, chk1, tst1_z))
		return tst1;
	arc4(msg2, msg2_z);
	if (ret < 0)
		return msg2;
	varg = (char **)calloc(argc + 10, sizeof(char *));
	if (!varg)
		return 0;
	if (ret) {
		arc4(rlax, rlax_z);
		if (!rlax[0] && key_with_file(shll))
			return shll;
		arc4(opts, opts_z);
		arc4(text, text_z);
		arc4(tst2, tst2_z);
		 key(tst2, tst2_z);
		arc4(chk2, chk2_z);
		if ((chk2_z != tst2_z) || memcmp(tst2, chk2, tst2_z))
			return tst2;
		/* Prepend hide_z spaces to script text to hide it. */
		scrpt = malloc(hide_z + text_z);
		if (!scrpt)
			return 0;
		memset(scrpt, (int) ' ', hide_z);
		memcpy(&scrpt[hide_z], text, text_z);
	} else {			/* Reexecute */
		if (*xecc) {
			scrpt = malloc(512);
			if (!scrpt)
				return 0;
			sprintf(scrpt, xecc, me);
		} else {
			scrpt = me;
		}
	}
	j = 0;
#if BUSYBOXON
	varg[j++] = "busybox";
	varg[j++] = "sh";
#else
	varg[j++] = argv[0];		/* My own name at execution */
#endif
	if (ret && *opts)
		varg[j++] = opts;	/* Options on 1st line of code */
	if (*inlo)
		varg[j++] = inlo;	/* Option introducing inline code */
	varg[j++] = scrpt;		/* The script itself */
	if (*lsto)
		varg[j++] = lsto;	/* Option meaning last option */
	i = (ret > 1) ? ret : 0;	/* Args numbering correction */
	while (i < argc)
		varg[j++] = argv[i++];	/* Main run-time arguments */
	varg[j] = 0;			/* NULL terminated array */
#if DEBUGEXEC
	debugexec(shll, j, varg);
#endif
	execvp(shll, varg);
	return shll;
}

int main(int argc, char ** argv)
{
#if DEBUGEXEC
	debugexec("main", argc, argv);
#endif
#if !TRACEABLE
	untraceable(argv[0]);
#endif
	argv[1] = xsh(argc, argv);
	fprintf(stderr, "%s%s%s: %s\n", argv[0],
		errno ? ": " : "",
		errno ? strerror(errno) : "",
		argv[1] ? argv[1] : "<null>"
	);
	return 1;
}
