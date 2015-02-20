#ifndef _o_Notify_h
#define _o_Notify_h


#define  O_NFY_ALWAYS    0
#define  O_NFY_FATAL     1
#define  O_NFY_WARN      2
#define  O_NFY_NOTICE    3
#define  O_NFY_INFO      4
#define  O_NFY_DEBUG     5
#define  O_NFY_DEBUG_01  6
#define  O_NFY_DEBUG_02  7
#define  O_NFY_DEBUG_03  8
#define  O_NFY_DEBUG_04  9 
#define  O_NFY_DEBUG_05  10
#define  O_NFY_NEVER     100

#define O_NFY_MSG_LENGTH  512

#define O_NFY_ENV  "NFY_LEVEL"


#ifdef __cplusplus
extern "C" {
#endif
  
  /** \file o_Notify.h
      print message. Usage is like printf() plus the print level as first parameter.
      
      The environment variable NFY_LEVEL can be set to override the value
      specified in o_NotifyLevel (i.e. for tcsh: "setenv NFY_LEVEL warn" for level O_NFY_WARN,
      the value string of the envirinment isn't case sensitiv, so warn is eq. to WARN). 
      If the level is set via NFY_LEVEL it can not be changed by an application. The seeting 
      with environment has higher priority than the level set by the application. Default level is
      O_NFY_INFO.  
  */
  
  
  /*!
     print message. 
     Usage is like printf() plus the print level as first parameter.
   */
  extern void o_Notify (int level, const char *format, ...);
  
  /*!   set the print level and can be one of the following. The print levels are sorted in increasing order. All messages are printed if they are equal or below current print level.

    \htmlonly
    <table border="1" cellspacing="4" cellpadding="9">
    <tr align="center"> <th> level </th> <th> description</th> <th> value </th> </tr> 
    <tr align="center"><td> <code>O_NFY_ALWAYS</code></td>   <td>always print to stdout</td><td><code>ALWAYS</code></td></tr> 
    <tr align="center"><td> <code>O_NFY_FATAL</code></td>    <td>for fatal error prints to stderr</td><td><code>FATAL</code></td></tr>
    <tr align="center"><td> <code>O_NFY_WARN</code></td>     <td>warning printing to stderr</td><td><code>WARN</code></td></tr>  
    <tr align="center"><td> <code>O_NFY_NOTICE</code></td>   <td> print notes to stdout</td><td><code>NOTICE</code></td></tr>
    <tr align="center"><td> <code>O_NFY_INFO</code></td>     <td> print info to stdout </td><td><code>INFO</code></td></tr>
    <tr align="center"><td> <code>O_NFY_DEBUG</code></td>    <td> print debug to stdout (general level)</td><td><code>DEBUG</code></td></tr>
    <tr align="center"><td> <code>O_NFY_DEBUG_01</code></td> <td> print debug to stdout (lowest level)</td><td><code>DEBUG_01</code></td></tr>
    <tr align="center"><td> <code>O_NFY_DEBUG_02</code></td> <td> print debug to stdout</td><td><code>DEBUG_02</code></td></tr>
    <tr align="center"><td> <code>O_NFY_DEBUG_03</code></td>  <td> print debug to stdout</td><td><code>DEBUG_03</code></td></tr>
    <tr align="center"><td> <code>O_NFY_DEBUG_04</code></td> <td> print debug to stdout</td><td><code>DEBUG_04</code></td></tr>
    <tr align="center"><td> <code>O_NFY_DEBUG_05</code></td> <td> print debug to stdout (highest level) </td><td><code>DEBUG_5</code></td></tr>   
    <tr align="center"><td> <code>O_NFY_NEVER</code></td> <td> don't print </td><td><code>NEVER</code></td></tr>
    </table>
    <p>
    \endhtmlonly
   */
  extern void o_NotifyLevel (int level);
  
  
  /*!   returns the print level
   */
  extern int  o_GetNotifyLevel(void);
  
  /*! trues the encryption on or off
  */
  extern void o_SetNotifyEncryption(bool on);
  
  /*! sets the encryption keyword
  */
  extern void o_SetNotifyEncryptionKeyword(char *key);
  
#ifdef __cplusplus
}
#endif


#endif /* _o_Notify */
