#ifndef __AUXLIB_H__
#define __AUXLIB_H__

#include <stdarg.h>

/* set_flags -
 *    Takes a string argument, and sets a flag for each char in the
 *    string.  As a special case, '@', sets all flags.
 * get_flag -
 *    Used by the DEBUGF macro to check to see if a flag has been set.
 *    Not supposed to be called directly by user code.
 */

void debug_set_flags(const char *optflags);
int debug_get_flag(char flag);

void debug_where(char flag, const char *file, int line,
                 const char *pretty_function, const char *msg, ...);
void debug_where_short(char flag, const char *file, int line, const char *msg,
                       ...);

/* DEBUGL "long debug":
 * Macro which expands into debug code.  First argument is a
 * debug flag char, second argument is output code that can
 * be sandwiched between <<.  Beware of operator precedence.
 * For example:
 *      DEBUGL ('u', "foo = " << foo);
 * will print two words and a newline if flag 'u' is  on.
 * Traces are preceded by filename, line number, and function.
 *
 * DEBUGS "short debug":
 * DEBUGL, but shortens debug statments by not printing the
 * long __PRETTY_FUNCTION__ macro
 *
 * DEBUGE "expression debug":
 * Macro which expands into debug code. Prints identical information to
 * DEBUGS, but without a message. Instead, one or more lines of code, which
 * are enclosed in the braces, are executed if the appropriate flag is set.
 *
 * If you need to do some work to get the information in a format that's
 * useful for debug output, you should put it in a DEBUGE with the same
 * flag as the message so it doesn't show up when debugging is disabled.
 */

#ifdef NDEBUG
#define DEBUGL(FLAG, PRINTF_ARGS...) ;
#define DEBUGE(FLAG, PRINTF_ARGS...) ;
#define DEBUGS(FLAG, PRINTF_ARGS...) ;
#else
#define DEBUGL(FLAG, PRINTF_ARGS...)                                           \
  {                                                                            \
    if (debug_get_flag(FLAG)) {                                                \
      debug_where(FLAG, __FILE__, __LINE__, __PRETTY_FUNCTION__, PRINTF_ARGS); \
    }                                                                          \
  }
#define DEBUGE(FLAG, STMT)                                            \
  {                                                                   \
    if (debug_get_flag(FLAG)) {                                       \
      debug_where(FLAG, __FILE__, __LINE__, __PRETTY_FUNCTION__, ""); \
      STMT;                                                           \
    }                                                                 \
  }
#define DEBUGS(FLAG, PRINTF_ARGS...)                            \
  {                                                             \
    if (debug_get_flag(FLAG)) {                                 \
      debug_where_short(FLAG, __FILE__, __LINE__, PRINTF_ARGS); \
    }                                                           \
  }
#endif
#endif
