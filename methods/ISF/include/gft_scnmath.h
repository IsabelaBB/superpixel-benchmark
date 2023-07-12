
#ifndef _GFT_SCNMATH_H_
#define _GFT_SCNMATH_H_

#include "gft_common.h"
#include "gft_scene.h"

namespace gft{
  namespace Scene32{

    Scene32 *Sub(Scene32 *scn1, Scene32 *scn2);
    /**
     * Inplace version.
     */
    void     Subinplace(Scene32 *scn1, Scene32 *scn2);

    Scene32 *Add(Scene32 *scn1, Scene32 *scn2);
    Scene32 *Add(Scene32 *scn, int value);
    /**
     * Inplace version.
     */
    void     Addinplace(Scene32 *scn1, Scene32 *scn2);

    Scene32 *Mult(Scene32 *scn1, Scene32 *scn2);

    Scene32 *Or(Scene32 *scn1, Scene32 *scn2);
    /**
     * Inplace version.
     */
    void     Orinplace(Scene32 *scn1, Scene32 *scn2);
    Scene32 *And(Scene32 *scn1, Scene32 *scn2);
    Scene32 *XOr(Scene32 *scn1, Scene32 *scn2);

    Scene32 *Complement(Scene32 *scn);
    Scene32 *Abs(Scene32 *scn);
    /**
     * Inplace version.
     */
    void     Negateinplace(Scene32 *scn);

  } //end Scene32 namespace


  namespace Scene16{

    Scene16 *Sub(Scene16 *scn1, Scene16 *scn2);
    /**
     * Inplace version.
     */
    void     Subinplace(Scene16 *scn1, Scene16 *scn2);

    Scene16 *Add(Scene16 *scn1, Scene16 *scn2);
    Scene16 *Add(Scene16 *scn, ushort value);
    /**
     * Inplace version.
     */
    void     Addinplace(Scene16 *scn1, Scene16 *scn2);

    Scene16 *Mult(Scene16 *scn1, Scene16 *scn2);

    Scene16 *Or(Scene16 *scn1, Scene16 *scn2);
    /**
     * Inplace version.
     */
    void     Orinplace(Scene16 *scn1, Scene16 *scn2);
    Scene16 *And(Scene16 *scn1, Scene16 *scn2);
    Scene16 *XOr(Scene16 *scn1, Scene16 *scn2);

    Scene16 *Complement(Scene16 *scn);

  } //end Scene16 namespace


  namespace Scene8{

    Scene8 *Sub(Scene8 *scn1, Scene8 *scn2);
    /**
     * Inplace version.
     */
    void    Subinplace(Scene8 *scn1, Scene8 *scn2);

    Scene8 *Add(Scene8 *scn1, Scene8 *scn2);
    Scene8 *Add(Scene8 *scn, uchar value);
    /**
     * Inplace version.
     */
    void    Addinplace(Scene8 *scn1, Scene8 *scn2);

    Scene8 *Mult(Scene8 *scn1, Scene8 *scn2);

    Scene8 *Or(Scene8 *scn1, Scene8 *scn2);
    /**
     * Inplace version.
     */
    void    Orinplace(Scene8 *scn1, Scene8 *scn2);
    Scene8 *And(Scene8 *scn1, Scene8 *scn2);
    Scene8 *XOr(Scene8 *scn1, Scene8 *scn2);

    Scene8 *Complement(Scene8 *scn);

  } //end Scene8 namespace


  namespace Scene{

    Scene *Sub(Scene *scn1, Scene *scn2);
    /**
     * Inplace version.
     */
    void   Subinplace(Scene *scn1, Scene *scn2);

    Scene *Add(Scene *scn1, Scene *scn2);
    Scene *Add(Scene *scn, int value);
    /**
     * Inplace version.
     */
    void   Addinplace(Scene *scn1, Scene *scn2);

    Scene *Mult(Scene *scn1, Scene *scn2);

    Scene *Or(Scene *scn1, Scene *scn2);
    /**
     * Inplace version.
     */
    void   Orinplace(Scene *scn1, Scene *scn2);
    Scene *And(Scene *scn1, Scene *scn2);
    Scene *XOr(Scene *scn1, Scene *scn2);

    Scene *Complement(Scene *scn);

  } //end Scene namespace


} //end gft namespace

#endif

