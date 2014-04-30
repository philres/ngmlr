#ifndef __ICONFIG_H__
#define __ICONFIG_H__

static char const * const MATCH_BONUS = 		"match_bonus";
static char const * const MATCH_BONUS_TT = 		"match_bonus_tt";
static char const * const MATCH_BONUS_TC = 		"match_bonus_tc";
static char const * const MISMATCH_PENALTY = 	"mismatch_penalty";
static char const * const GAP_REF_PENALTY = 	"gap_ref_penalty";
static char const * const GAP_READ_PENALTY = 	"gap_read_penalty";
static char const * const GAP_EXTEND_PENALTY = 	"gap_extend_penalty";

static char const * const MODE =			 	"mode";

static char const * const MIN_MQ =			 	"min_mq";

static char const * const LOCAL =			 	"local";
static char const * const ENDTOEND =		 	"end_to_end";

static char const * const KEEPTAGS =		 	"keep_tags";

static char const * const SKIP_MATE_CHECK =		"skip_mate_check";

static char const * const RG_ID =		 	    "rg_id";
static char const * const RG_CN =		 	    "rg_cn";
static char const * const RG_DS =		 	    "rg_ds";
static char const * const RG_DT =		 	    "rg_dt";
static char const * const RG_FO =		 	    "rg_fo";
static char const * const RG_KS =		 	    "rg_ks";
static char const * const RG_LB =		 	    "rg_lb";
static char const * const RG_PG =		 	    "rg_pg";
static char const * const RG_PI =		 	    "rg_pi";
static char const * const RG_PL =		 	    "rg_pl";
static char const * const RG_PU =		 	    "rg_pu";
static char const * const RG_SM =		 	    "rg_sm";

static char const * const MAX_C_INDELS =		 	    "max_consec_indels";


#ifdef DEBUGLOG
static char const * const LOG =			 	"log";
static char const * const LOG_LVL =			 	"log_lvl";
#endif


class IConfig
{
public:

	virtual char const * GetString(char const * const name) const = 0;
	virtual int GetInt(char const * const name) const = 0;
	virtual int GetInt(char const * const name, int min, int max) const = 0;
	virtual int GetParameter(char const * const name) const = 0;
	virtual float GetFloat(char const * const name) const = 0;
	virtual float GetFloat(char const * const name, float min, float max) const = 0;

	virtual int GetIntArray(char const * const name, int * pData, int len) const = 0;
	virtual int GetFloatArray(char const * const name, float * pData, int len) const = 0;
	virtual int GetDoubleArray(char const * const name, double * pData, int len) const = 0;

	virtual bool Exists(char const * const name) const = 0;
	virtual bool HasArray(char const * const name) const = 0;

	virtual ~IConfig() {};
};

typedef void (*pfSetConfig)(IConfig const *);

extern IConfig* _config;
#define Config (*_config)

#endif
