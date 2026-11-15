"""
Turn solver output (2D profiles through time) into a mesh for the 3D viewer.

Mesh columns  = water particles (material labels), so the grid follows the real flow.
Mesh rows     = time.  A wave meeting a straight reef at an angle is the same 2D event
                delayed along the reef, so time maps to distance along the reef.
After the lip lands potential flow ends; rows after touchdown are a modelled bore
(clearly not solver output) so the whitewater side of the picture is filled in.
"""
import numpy as np, sys, json, base64
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d

g = 9.81
name = sys.argv[1]
d = np.load(name + '.npz')
T = d['t']; n = d['n']; off = np.concatenate([[0], np.cumsum(n)])
F = [dict(t=T[i], **{k: d[k][off[i]:off[i + 1]] for k in 'xyauv'}) for i in range(len(T))]
bx, by = d['bx'], d['by']
tb = float(T[-1])
Tp, Ta, Tw = 14.0, 9.0, 5.0
ts, te = tb - Ta, tb + Tw
sm = lambda a, b, x: np.clip((x - a) / (b - a), 0, 1) ** 2 * (3 - 2 * np.clip((x - a) / (b - a), 0, 1))

# ---------- landmarks on the touchdown frame
f = F[-1]
x, y, a = f['x'], f['y'], f['a']
s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
ic = int(np.argmax(y))
j = ic
while x[j + 1] >= x[j]: j += 1
itip = j
while x[j + 1] <= x[j]: j += 1
iback = j
while x[j] < x[itip]: j += 1
itd = j
front = (x > x[ic]) & (x < x[ic] + 20)
ytr = y[front].min()
Hb = y[ic] - ytr
xtd = x[itip]
print(f"{name}: tb={tb:.2f} crest y={y[ic]:.2f} trough={ytr:.2f} H={Hb:.2f} tip=({x[itip]:.2f},{y[itip]:.2f}) tube back x={x[iback]:.2f}")

# ---------- choose columns (labels) from the touchdown frame
A0x, A1x = 30.0, xtd + 60.0
sA = np.interp(x[ic] - 1.3 * Hb, x[:ic], s[:ic]); sB = s[itd] + 1.0 * Hb
ds_b, ds_far, grade = 0.11, 2.6, 0.10
sf = np.arange(s[np.searchsorted(x[:ic], A0x)], s[itd + np.searchsorted(x[itd:], A1x)], 0.005)
dist = np.maximum(0, np.maximum(sA - sf, sf - sB))
hh = np.minimum(ds_far, ds_b + grade * dist)
hh = np.minimum(hh, 0.035 + 0.12 * np.abs(sf - s[itip]))
nu = np.concatenate([[0], np.cumsum(0.005 / hh)])[:-1]
M = int(nu[-1]) + 1
scol = np.interp(np.arange(M), nu, sf)
k = np.argmin(np.abs(scol - s[itip])); scol += s[itip] - scol[k]          # a column exactly on the tip
acol = np.interp(scol, s, a)
ctip = k
ctd = int(np.argmin(np.abs(scol - s[itd]))); cc = int(np.argmin(np.abs(scol - s[ic])))
print("columns", M, "tip col", ctip, "td col", ctd)

def sample(fr, lab):
    """positions / velocities of the labelled particles in one solver frame"""
    idx = np.arange(len(fr['a']))
    fi = np.interp(lab, fr['a'], idx)
    out = {}
    da = CubicSpline(idx, fr['a'])(fi, 1)
    for kk in 'xyuv':
        sp = CubicSpline(idx, fr[kk])
        out[kk] = sp(fi)
        if kk in 'xy': out['d' + kk] = sp(fi, 1) / da
    return out

def sample_t(t, lab):
    i = np.searchsorted(T, t)
    if i <= 0: return sample(F[0], lab)
    if i >= len(T): return sample(F[-1], lab)
    w = (t - T[i - 1]) / (T[i] - T[i - 1])
    p, q = sample(F[i - 1], lab), sample(F[i], lab)
    return {kk: (1 - w) * p[kk] + w * q[kk] for kk in p}

# ---------- rows before touchdown: straight from the solver
rows_t = [tb]
while rows_t[-1] > ts + 1e-6:
    t = rows_t[-1]
    dt = float(np.clip(0.03 + 0.10 * (tb - t - 0.35), 0.02, 0.1))
    i = np.argmin(np.abs(T - (t - dt)))
    if T[i] >= t - 1e-6: i -= 1
    rows_t.append(max(float(T[i]), ts) if T[i] > ts else ts)
rows_t = rows_t[::-1]
pre = [sample_t(t, acol) for t in rows_t]
npre = len(pre)
fin = pre[-1]

# ---------- rows after touchdown: modelled bore
post_d = np.concatenate([[0.02, 0.04, 0.07, 0.10, 0.14, 0.18, 0.23, 0.28, 0.34, 0.4], np.arange(0.5, 1.6, 0.1),
                         np.arange(1.75, Tw + 1e-6, 0.25)])
x0, y0 = fin['x'].copy(), fin['y'].copy()
utip = float(fin['u'][ctip])
# wave speed from the motion of the back wall of the tube over the last half second
pb = sample_t(tb - 0.5, acol)
cw = float((fin['x'][cc] - pb['x'][cc]) / 0.5)
cinf = 0.8 * cw
xf = lambda D: xtd + cinf * D + (utip - cinf) * 0.7 * (1 - np.exp(-D / 0.7))
Dg = np.linspace(0, 40, 4001); xfg = xf(Dg)
Dpass = lambda X: np.interp(X, xfg, Dg)
hb = lambda D: Hb * (0.50 * np.exp(-D / 2.4) * (1 - np.exp(-D / 0.12)) + 0.42 * (D / 0.32) * np.exp(1 - D / 0.32))
print(f"utip={utip:.2f} cw={cw:.2f}")
Ax, Ay = x0[:ctip + 1], y0[:ctip + 1]                 # back + lip top: a height field
Cx, Cy = x0[ctd:], y0[ctd:]                           # water ahead of the landing point
wA = sm(-75., -16., Ax - xtd)
jj = np.arange(M)

def target(D):
    shift = xf(D) - xtd
    X = np.empty(M); Y = np.empty(M)
    XA = Ax + wA * shift
    lvl = np.interp(xf(D), Cx, Cy) * np.exp(-D / 2.5)
    E = lvl + (np.interp(XA - shift, Ax, Ay) - ytr) * np.exp(-D / 2.0)
    # far behind the wave the envelope is not shifted: blend to the un-shifted surface
    E = (1 - wA) * Ay + wA * E
    plat = lvl + hb(D)
    kk = 0.5
    smx = kk * np.log(np.exp((E - lvl) / kk) + np.exp((plat - lvl) / kk) - 1 + 1e-12) + lvl
    # the plateau only exists near the front; fade it out over ~3 bore heights behind
    xi = xf(D) - XA
    near = np.exp(-np.maximum(xi, 0) / (2.0 + 3.0 * hb(D)))
    top = E + wA * near * (smx - E)
    w = 0.25 + 1.5 * hb(D) * min(1.0, D / 0.45)
    q = np.clip(xi / w, 0, 1)
    nose = np.clip(1 - (1 - q) ** 2, 0, 1) ** 0.6
    X[:ctip + 1] = XA; Y[:ctip + 1] = lvl + (top - lvl) * nose
    X[ctip + 1:ctd] = xf(D); Y[ctip + 1:ctd] = lvl
    XC = np.maximum(Cx, xf(D) + 0.003 * (jj[ctd:] - ctd))
    X[ctd:] = XC; Y[ctd:] = np.interp(XC, Cx, Cy) * np.exp(-D / 2.5)
    return X, Y

post = []
for D in post_d:
    X, Y = target(D)
    if D < 0.36:
        wgt = sm(0, 0.36, D)
        bxp = x0 + fin['u'] * D; byp = y0 + fin['v'] * D - 0.5 * g * D * D * sm(0, 1, (y0 - ytr) / 1.0) * 0
        floor = np.interp(bxp, Cx, Cy, left=-1e9)
        byp = np.maximum(byp, floor - 0.02)
        X = (1 - wgt) * bxp + wgt * X; Y = (1 - wgt) * byp + wgt * Y
    wend = 1 - sm(Tw - 2.0, Tw, D)
    nxt = sample_t(tb + D - Tp, acol)                 # the next wave of the set, arriving offshore
    X = acol + wend * (X - acol) + (nxt['x'] - acol)
    Y = wend * Y + nxt['y']
    post.append(dict(x=X, y=Y))
rows_t = list(rows_t) + [tb + D for D in post_d]
R = len(rows_t)
PX = np.array([p['x'] for p in pre] + [p['x'] for p in post])
PY = np.array([p['y'] for p in pre] + [p['y'] for p in post])
tau = np.array(rows_t)
print("rows", R, "pre", npre, "verts", R * M)
assert np.abs(PX[-1] - PX[0]).max() < 1e-6 and np.abs(PY[-1] - PY[0]).max() < 1e-6

# velocities (d/dtau at fixed column) and tangents (d/dcolumn)
VU = np.zeros_like(PX); VV = np.zeros_like(PX)
for r in range(R):
    if r < npre:
        VU[r], VV[r] = pre[r]['u'], pre[r]['v']
    else:
        r0, r1 = r - 1, min(r + 1, R - 1)
        if r == R - 1:
            VU[r], VV[r] = VU[0], VV[0]; continue
        VU[r] = (PX[r1] - PX[r0]) / (tau[r1] - tau[r0]); VV[r] = (PY[r1] - PY[r0]) / (tau[r1] - tau[r0])
TX = np.gradient(PX, axis=1); TY = np.gradient(PY, axis=1)
for r in range(npre):
    TX[r], TY[r] = pre[r]['dx'], pre[r]['dy']
da_col = np.gradient(acol)
for r in range(npre, R):
    TX[r] /= da_col; TY[r] /= da_col
STR = np.hypot(TX, TY)                                # stretch: metres of surface per metre of label
STR = np.maximum(STR, 1e-4)
TX /= STR; TY /= STR
# collapsed columns after touchdown have no direction of their own: borrow a level tangent
bad = (np.hypot(np.gradient(PX, axis=1), np.gradient(PY, axis=1)) < 1e-5)
TX[bad] = 1; TY[bad] = 0

# ---------- lighting helpers computed in the profile plane
SUN_EL, SUN_AZ = np.radians(40.0), np.radians(28.0)
sun3 = np.array([-np.cos(SUN_EL) * np.cos(SUN_AZ), np.sin(SUN_EL), np.cos(SUN_EL) * np.sin(SUN_AZ)])
sang = np.arctan2(sun3[1], sun3[0]); s2 = np.array([np.cos(sang), np.sin(sang)])
K = 60
angs = (np.arange(K) + 0.5) * np.pi / K
dirs = np.stack([np.cos(angs), np.sin(angs)], 1)
TH1 = np.zeros_like(PX); TH2 = np.full_like(PX, np.pi); THK = np.full_like(PX, 30.0); SHD = np.zeros_like(PX)

def cast(P, sel, D):
    """first hit distance of rays P[sel] + t*D[k] against polyline P.  returns (nsel, nk)"""
    A = P[:-1]; B = P[1:]; Eg = B - A
    O = P[sel]
    res = np.full((len(sel), len(D)), np.inf)
    for ki, dv in enumerate(D):
        den = dv[0] * Eg[:, 1] - dv[1] * Eg[:, 0]                       # (m,)
        AO = A[None, :, :] - O[:, None, :]                              # (n,m,2)
        with np.errstate(divide='ignore', invalid='ignore'):
            tt = (AO[..., 0] * Eg[None, :, 1] - AO[..., 1] * Eg[None, :, 0]) / den[None, :]
            uu = (AO[..., 0] * dv[1] - AO[..., 1] * dv[0]) / den[None, :]
        ok = (tt > 0.02) & (uu >= 0) & (uu <= 1)
        seg = np.arange(len(A))[None, :]
        ok &= (np.abs(seg - sel[:, None]) > 1) & (np.abs(seg + 1 - sel[:, None]) > 1)
        tt = np.where(ok, tt, np.inf)
        res[:, ki] = tt.min(1)
    return res

LIT = np.ones_like(PX)
sel = np.arange(M)
for r in range(R - 1):
    P = np.stack([PX[r], PY[r]], 1)
    nrm = np.stack([-TY[r], TX[r]], 1)
    facing = nrm @ s2
    hs = cast(P, sel, s2[None, :])[:, 0]
    wf = sm(-0.20, 0.25, facing)                       # 0 = facing away from the sun, 1 = facing it
    th = np.where(np.isfinite(hs), np.minimum(hs, 24.0), 24.0)
    THK[r] = wf * 24.0 + (1 - wf) * th
    hd = cast(P, sel, np.array([[0.0, -1.0]]))[:, 0]           # thin sheet of water with air under it (the lip)
    thin = np.isfinite(hd) & (nrm[:, 1] > 0.05)
    THK[r] = np.where(thin, np.minimum(THK[r], 1.2 + 1.6 * hd), THK[r])
    lit_away = 0.40 + 0.60 * np.exp(-0.17 * th)
    lit_to = np.where(np.isfinite(hs), 0.32, 1.0)
    LIT[r] = wf * lit_to + (1 - wf) * lit_away
    SHD[r] = wf * np.where(np.isfinite(hs), 1.0, 0.0)
    steep = np.abs(TY[r]).max() > 0.5 or tau[r] > tb
    if not steep or tau[r] > tb + 0.6: continue
    hit = cast(P, sel, dirs)
    free = ~np.isfinite(hit)
    for c in sel:
        fr = free[c]
        if fr.all(): continue
        best = (0, 0); st = None
        for kk in range(K + 1):
            if kk < K and fr[kk]:
                if st is None: st = kk
            else:
                if st is not None and kk - st > best[1] - best[0]: best = (st, kk)
                st = None
        TH1[r, c] = best[0] * np.pi / K; TH2[r, c] = best[1] * np.pi / K
THK[-1] = THK[0]; LIT[-1] = LIT[0]; SHD[-1] = SHD[0]
for arr, sg in ((TH1, 1.2), (TH2, 1.2), (SHD, 1.5), (LIT, 1.5), (THK, 1.0)):
    arr[:] = gaussian_filter1d(arr, sg, axis=1)

# ---------- foam
FOAM = np.zeros_like(PX); TURB = np.zeros_like(PX)
arc = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x0), np.diff(y0)))])
for r in range(R):
    D = tau[r] - tb
    if D <= 0:
        # lip edge feathering as the jet runs out
        reach = 2.4 * sm(-0.55, 0.0, D) ** 1.6
        P = np.stack([PX[r], PY[r]], 1)
        ar = np.concatenate([[0], np.cumsum(np.hypot(np.diff(PX[r]), np.diff(PY[r])))])
        dtip = np.abs(ar - ar[ctip])
        FOAM[r] = np.where(reach > 0.02, np.clip(1 - dtip / max(reach, 1e-3), 0, 1) ** 1.2, 0)
        # foam left by the previous wave of the set
        ap = D + Tp - Dpass(PX[r])
        FOAM[r] += np.where((PX[r] > xtd - 1) & (ap > 0), 0.9 * np.exp(-ap / 3.0), 0)
    else:
        X = PX[r]; xi = xf(D) - X
        age = np.where(X >= xtd, D - Dpass(np.clip(X, xtd, None)), D)
        behind = np.where(X < xtd, np.exp(-(xtd - X) / (2.0 + 2.5 * min(D, 1.5))), 1.0)
        on = sm(-0.05, 0.12, xi)
        FOAM[r] = on * behind * np.exp(-np.maximum(age, 0) / 4.0)
        TURB[r] = on * behind * (0.10 + 0.22 * hb(D) / Hb * 2.2) * np.exp(-np.maximum(age, 0) / 1.3) * (1 - sm(Tw - 2.0, Tw, D)) * sm(0.0, 0.25, D)
        ap = D + Tp - Dpass(X)
        FOAM[r] += np.where((X > xf(D) + 0.2) & (X > xtd), 0.9 * np.exp(-ap / 3.0), 0)
FOAM = np.clip(FOAM, 0, 1)
DEP = PY - np.interp(PX, bx, by)

# ---------- pack
def q16(arr, lo, hi):
    return np.round((np.clip(arr, lo, hi) - lo) / (hi - lo) * 65535).astype('<u2')
umax = float(np.abs(np.stack([VU, VV])).max()) + 0.01
smax = 40.0
bufs = [PX.astype('<f4'), PY.astype('<f4'),
        q16(TX, -1, 1), q16(TY, -1, 1), q16(VU, -umax, umax), q16(VV, -umax, umax),
        q16(np.log(STR), np.log(1e-2), np.log(smax)), q16(TH1, 0, np.pi), q16(TH2, 0, np.pi),
        q16(THK, 0, 30), q16(SHD, 0, 1), q16(LIT, 0, 1), q16(FOAM, 0, 1), q16(TURB, 0, 2), q16(DEP, 0, 16)]
raw = b''.join(b.tobytes() for b in bufs)
trk = lambda arr: [round(float(v), 3) for v in arr]
meta = dict(name=name, M=M, R=R, npre=npre, tb=tb, Tp=Tp, Tw=Tw, umax=umax, smax=smax,
            tau=trk(tau), acol=trk(acol), ctip=int(ctip), ctd=int(ctd), cc=int(cc),
            xtd=float(xtd), ytr=float(ytr), Hb=float(Hb), utip=utip, cw=cw, cinf=cinf,
            vtip=float(fin['v'][ctip]),
            tipx=trk(PX[:npre, ctip]), tipy=trk(PY[:npre, ctip]), tipu=trk(VU[:npre, ctip]), tipv=trk(VV[:npre, ctip]),
            crx=trk([PX[r, np.argmax(PY[r, :ctd])] for r in range(R)]), cry=trk([PY[r, :ctd].max() for r in range(R)]),
            sun=trk(sun3), bx=trk(bx[::6]), by=trk(by[::6]))
open(f'../data_{name}.js', 'w').write(
    f"window.WAVE_DATA=window.WAVE_DATA||{{}};window.WAVE_DATA[{json.dumps(name)}]="
    + json.dumps(dict(meta=meta, b64=base64.b64encode(raw).decode())) + ";")
print("bytes", len(raw), "umax", umax, "stretch range", STR.min(), STR.max())
np.savez(f'{name}_mesh.npz', LIT=LIT, PX=PX, PY=PY, tau=tau, FOAM=FOAM, TURB=TURB, TH1=TH1, TH2=TH2, THK=THK, SHD=SHD)
