#!/usr/bin/env python3
"""Evaluation wrapper for log-PDF spline (3D: dr, phiE, dt).

Needed for how this spline must be consumed.  The .fits stores
log(PDF) and is unit-normalised over dt per (dr, phiE) slice, so a correct lookup must do the following:
  1. clamp (dr, phiE) into spline.extents - a coordinate past the table edge uses the
     nearest valid slice
  2. FLOOR out-of-support dt to LOG_FLOOR: photospline returns 0
     in log-space outside its knots, and exp(0)=1 would otherwise read as the MOST
     probable time. Flooring to LOG_FLOOR puts
     it at the bottom of the distribution instead.
  3. exp() back to density.
  4. (evalPdf only) renormalise to unit dt-area

Both validation and reco scripts should import this, so the validated curves are exactly the curve
reco evaluates.
"""

import numpy as np

LOG_FLOOR = -30.0


def evalLogPdf(spline, dr, phiE, dt):
    """PDF (density) at fixed scalar (dr, phiE) over a dt array-like.

    Clamps dr/phiE into extent, floors out-of-dt-support to LOG_FLOOR, then exp()."""
    # unpack all 3 ranges
    (drLo, drHi), (peLo, peHi), (dtLo, dtHi) = spline.extents

    # clamp dr and phiE
    drC = float(np.clip(dr, drLo, drHi))
    peC = float(np.clip(phiE, peLo, peHi))

    # normalize dt into a 1D float array.
    # Just as a safeguard because evaluate_simple takes an array (I think)
    dt = np.atleast_1d(np.asarray(dt, dtype=float))

    # clamp dt
    dtC = np.clip(dt, dtLo, dtHi)  # keep evaluate_simple supported

    # query the spline at the clamped coordinate
    # .squeeze drops extra dims
    logp = spline.evaluate_simple([[drC], [peC], dtC]).squeeze()

    # for each dt point, ask if it was originally inside the fitted window
    # if yes, keep the spline's value
    # if no, floor it to -30
    logp = np.where((dt >= dtLo) & (dt <= dtHi), logp, LOG_FLOOR)

    # return linear density (convert from log space)
    return np.exp(logp)


def evalPdf(spline, dr, phiE, dt):
    """Unit-normalised PDF: evalLogPdf divided by its own dt-integral.

    dt must span the spline's dt support for the normalisation to be meaningful
    (e.g. np.linspace(dtLo, dtHi, N))."""
    dt = np.atleast_1d(np.asarray(dt, dtype=float))
    pdf = evalLogPdf(spline, dr, phiE, dt)
    area = np.trapz(pdf, dt)
    return pdf / area if area > 0 else pdf
