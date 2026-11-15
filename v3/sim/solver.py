"""
2D fully-nonlinear potential-flow solver for an overturning wave over a reef.

Method: mixed Eulerian-Lagrangian boundary integral (dipole / vortex-sheet form,
after Longuet-Higgins & Cokelet 1976, Baker-Meiron-Orszag 1982).
 - fluid is inviscid, incompressible, irrotational: laplace(phi) = 0
 - free surface markers move with the water: Dz/Dt = grad(phi)
 - Bernoulli on the surface: Dphi/Dt = |u|^2/2 - g y
 - periodic in x (period L); the reef is put in with a conformal map z = f(zeta)
   from a flat-bottomed strip, so the bottom boundary condition is exact (images).
"""
import numpy as np, sys, json, time
from scipy.interpolate import CubicSpline
from scipy.ndimage import maximum_filter1d, gaussian_filter1d

g = 9.81


class Reef:
    """z = f(zeta) = zeta + sum b_n exp(i n k (zeta + i h)); flat bottom Im zeta = -h."""

    def __init__(self, L, depth_fn, M=48, smooth=3.0):
        self.L = L
        self.k = 2 * np.pi / L
        n = 2048
        xi = np.arange(n) * L / n
        d = depth_fn(xi)
        self.h = d.mean()
        Bt = (self.h - d)
        beta = np.fft.fft(Bt) / n
        nn = np.arange(1, M + 1)
        self.n = nn
        self.b = 2j * beta[1:M + 1] * np.exp(-(nn * self.k * smooth) ** 2 / 2)

    def f(self, zeta):
        E = np.exp(1j * self.k * np.multiply.outer(zeta + 1j * self.h, self.n))
        return zeta + E @ self.b

    def fp(self, zeta):
        E = np.exp(1j * self.k * np.multiply.outer(zeta + 1j * self.h, self.n))
        return 1 + E @ (1j * self.n * self.k * self.b)

    def inv(self, z):
        zeta = np.array(z, dtype=complex)
        for _ in range(60):
            dz = self.f(zeta) - z
            zeta = zeta - dz / self.fp(zeta)
            if np.abs(dz).max() < 1e-12:
                break
        return zeta

    def bottom(self, n=600):
        xi = np.arange(n + 1) * self.L / n
        return self.f(xi - 1j * self.h)


def dfft(p, order=1):
    N = len(p)
    kk = np.fft.fftfreq(N, 1.0 / N)
    P = np.fft.fft(p)
    if N % 2 == 0 and order % 2 == 1:
        P[N // 2] = 0
    return np.fft.ifft(P * (1j * kk) ** order)


def lowpass(p, strength=36.0, power=10):
    N = len(p)
    kk = np.abs(np.fft.fftfreq(N, 1.0 / N)) / (N / 2)
    return np.fft.ifft(np.fft.fft(p) * np.exp(-strength * kk ** power))


class Solver:
    def __init__(self, reef):
        self.reef = reef
        self.L = reef.L
        self.h = reef.h

    def derivs(self, zeta, phi):
        L, h = self.L, self.h
        N = len(zeta)
        de = 2 * np.pi / N
        e = np.arange(N) * de
        p = zeta - L * e / (2 * np.pi)
        ze = L / (2 * np.pi) + dfft(p, 1)
        zee = dfft(p, 2)
        D = (np.pi / L) * (zeta[None, :] - zeta[:, None])
        np.fill_diagonal(D, 1.0)
        C1 = 1.0 / np.tan(D)
        np.fill_diagonal(C1, 0.0)
        C2 = 1.0 / np.tan((np.pi / L) * (np.conj(zeta)[None, :] - 2j * h - zeta[:, None]))
        fac = de / (2j * L)
        K = (fac * (ze[None, :] * C1 - np.conj(ze)[None, :] * C2)).real
        K[np.diag_indices(N)] += (zee / ze).imag / (4 * np.pi) * de
        A = K
        A[np.diag_indices(N)] -= 0.5
        mu = np.linalg.solve(A, phi)
        mue = dfft(mu, 1).real
        muee = dfft(mu, 2).real
        c = mue / ze
        S1 = C1 @ mue - c * (C1 @ ze) + (L / np.pi) * (muee - c * zee) / ze
        S2 = C2 @ mue
        Wz = -0.5 * c + fac * (S1 - S2)
        fp = self.reef.fp(zeta)
        z = self.reef.f(zeta)
        q2 = np.abs(Wz) ** 2 / np.abs(fp) ** 2
        dzeta = np.conj(Wz) / np.abs(fp) ** 2
        dphi = 0.5 * q2 - g * z.imag
        return dzeta, dphi, Wz, ze, fp, z

    def step(self, zeta, phi, dt):
        k1z, k1p, *_ = self.derivs(zeta, phi)
        k2z, k2p, *_ = self.derivs(zeta + 0.5 * dt * k1z, phi + 0.5 * dt * k1p)
        k3z, k3p, *_ = self.derivs(zeta + 0.5 * dt * k2z, phi + 0.5 * dt * k2p)
        k4z, k4p, *_ = self.derivs(zeta + dt * k3z, phi + dt * k3p)
        return (zeta + dt / 6 * (k1z + 2 * k2z + 2 * k3z + k4z),
                phi + dt / 6 * (k1p + 2 * k2p + 2 * k3p + k4p))

    def energy(self, zeta, phi):
        dzeta, dphi, Wz, ze, fp, z = self.derivs(zeta, phi)
        N = len(zeta)
        de = 2 * np.pi / N
        zpe = ze * fp                      # dz/de physical
        u = np.conj(Wz / fp)               # u + i v
        # phi_n ds with outward normal (up): n ds = i * dz rotated -> (-dy, dx)
        phin_ds = (-(u.real) * zpe.imag + u.imag * zpe.real) * de
        KE = 0.5 * np.sum(phi * phin_ds)
        PE = 0.5 * g * np.sum(z.imag ** 2 * zpe.real) * de
        vol = np.sum(z.imag * zpe.real) * de
        return KE, PE, vol


def regrid(zeta, phi, lab, L, reef, ds_min=0.03, ds_max=1.2, cpr=0.14, Nmax=900, Nmin=256, grade=0.12, ds_face=0.08):
    """Redistribute markers: spacing ~ cpr * radius of curvature, clipped and smoothed."""
    N = len(zeta)
    e = np.arange(N + 1) * 2 * np.pi / N
    p = zeta - L * e[:-1] / (2 * np.pi)
    pp = np.append(p, p[0])
    spl_p = CubicSpline(e, pp, bc_type='periodic')
    spl_f = CubicSpline(e, np.append(phi, phi[0]), bc_type='periodic')
    spl_a = CubicSpline(e, np.append(lab, lab[0]), bc_type='periodic')
    # fine sampling
    M = 8 * N
    ef = np.arange(M) * 2 * np.pi / M
    zf = spl_p(ef) + L * ef / (2 * np.pi)
    ze = spl_p(ef, 1) + L / (2 * np.pi)
    zee = spl_p(ef, 2)
    sp = np.abs(ze)
    kap = np.abs((np.conj(ze) * zee).imag) / sp ** 3
    ds_f = sp * 2 * np.pi / M
    s = np.concatenate([[0], np.cumsum(ds_f)])
    Stot = s[-1]
    if not np.isfinite(Stot) or Stot > 2 * L:
        raise RuntimeError('surface blew up')
    dl = 0.01
    ns = int(Stot / dl)
    dl = Stot / ns
    su = (np.arange(ns) + 0.5) * dl
    kf = np.interp(su, 0.5 * (s[:-1] + s[1:]), kap, period=Stot)
    kf = gaussian_filter1d(kf, sigma=0.03 / dl, mode='wrap')
    hh = np.clip(cpr / (kf + 1e-9), ds_min, ds_max)
    # keep the whole breaking region (face, lip, tube) finely resolved once the face is steep
    zs = np.interp(su, s, np.append(zf.real, zf.real[0] + L)) + 1j * np.interp(su, s, np.append(zf.imag, zf.imag[0]))
    tx = np.gradient(zs.real, dl); ty = np.gradient(zs.imag, dl)
    if np.abs(ty).max() > 0.45:
        ic = np.argmax(zs.imag); xc = zs.real[ic]
        front = (zs.real > xc) & (zs.real < xc + 25)
        H = zs.imag[ic] - zs.imag[front].min()
        steep = np.abs(ty) > 0.45
        xa = zs.real[steep].min() - 0.6 * H; xb = zs.real[steep].max() + 0.8 * H
        hh = np.where((zs.real > xa) & (zs.real < xb), np.minimum(hh, ds_face), hh)
    # gradient-limited spacing: neighbouring cells differ by at most `grade` per unit length
    for _ in range(2):
        hh = np.minimum.accumulate(hh - grade * su) + grade * su
        hh = (np.minimum.accumulate((hh + grade * su)[::-1]))[::-1] - grade * su
        hh = np.roll(hh, ns // 2)
    nu = np.concatenate([[0], np.cumsum(dl / hh)])
    Nn = int(round(nu[-1]))
    Nn = min(max(Nn, Nmin), Nmax)
    Nn += Nn % 2
    # smooth the spacing in index space so that zeta(e) is a smooth function of e
    q = 8
    nuu = (np.arange(Nn * q) + 0.5) * nu[-1] / (Nn * q)
    hu = np.interp(nuu, 0.5 * (nu[:-1] + nu[1:]), hh, period=nu[-1])
    hu = gaussian_filter1d(hu, sigma=3.0 * q, mode='wrap')
    cs = np.concatenate([[0], np.cumsum(hu)])
    cs *= Stot / cs[-1]
    s_new = cs[::q][:Nn]
    e_new = np.interp(s_new, s, np.append(ef, 2 * np.pi))
    z_new = spl_p(e_new) + L * e_new / (2 * np.pi)
    return z_new, spl_f(e_new), spl_a(e_new)


def min_gap(z, sep=1.5):
    """smallest distance between surface points that are far apart along the surface"""
    N = len(z)
    ds = np.abs(np.diff(np.append(z, z[0] + 0)))
    ds[-1] = ds[-2]
    s = np.concatenate([[0], np.cumsum(ds)])[:-1]
    y = z.imag
    idx = np.where(y > y.min() + 0.25 * (y.max() - y.min()))[0]
    if len(idx) < 4:
        return 1e9
    zz = z[idx]
    ss = s[idx]
    Dm = np.abs(zz[:, None] - zz[None, :])
    far = np.abs(ss[:, None] - ss[None, :]) > sep
    if not far.any():
        return 1e9
    return Dm[far].min()


def run(cfg, out, verbose=True):
    L = cfg['L']
    dfn = cfg['depth']
    reef = Reef(L, dfn, M=cfg.get('M', 48), smooth=cfg.get('smooth', 3.0))
    sol = Solver(reef)
    N0 = cfg.get('N0', 320)
    x = np.arange(N0) * L / N0
    eta, u = cfg['init'](x)
    u = u - u.mean()
    phi = np.cumsum(u) * L / N0
    phi -= phi.mean()
    zeta = reef.inv(x + 1j * eta)
    lab = np.zeros(N0)          # label = x0 - L*e/2pi (periodic part), store periodic part
    lab0 = x.copy()
    labp = np.zeros(N0)         # periodic part of material label
    # the material label a = x0; store as a - L e/2pi so it stays periodic
    t = 0.0
    frames = []
    tout = 0.0
    dto = cfg.get('dto', 0.05)
    T = cfg['T']
    kw = dict(ds_min=cfg.get('ds_min', 0.03), ds_max=cfg.get('ds_max', 1.2), cpr=cfg.get('cpr', 0.14),
              Nmax=cfg.get('Nmax', 900), Nmin=cfg.get('Nmin', 256), grade=cfg.get('grade', 0.12), ds_face=cfg.get('ds_face', 0.08))
    step = 0
    E0 = None
    t0 = time.time()
    # label handling: keep full label a_j and interpolate (a - L e/(2 pi)) through regrids
    a = x.copy()
    gap_stop = cfg.get('gap', 0.12)
    import pickle, os
    ck = cfg.get('ckpt')
    if ck and cfg.get('resume') and os.path.exists(ck):
        t, zeta, phi, a, frames, E0, tout, step = pickle.load(open(ck, 'rb'))
        print('resumed at', t)
        ck_done = True
    else:
        ck_done = False
    while t < T:
        N = len(zeta)
        if ck and not ck_done and t >= cfg.get('ckpt_t', 1e9):
            pickle.dump((t, zeta, phi, a, frames, E0, tout, step), open(ck, 'wb'))
            ck_done = True
        if step % cfg.get('regrid_every', 4) == 0:
            try:
                zeta, phi, ap = regrid(zeta, phi, a - zeta.real, L, reef, **kw)
            except Exception as ex:
                print("regrid failed", ex); break
            a = ap + zeta.real
            N = len(zeta)
        dzeta, dphi, Wz, ze, fp, z = sol.derivs(zeta, phi)
        uu = np.conj(Wz / fp)
        dsz = np.abs(np.diff(np.append(z, z[0] + L)))
        rel = np.abs(np.diff(np.append(uu, uu[0])))
        dt = min(cfg.get('cfl_g', 0.35) * np.sqrt(dsz.min() / g),
                 cfg.get('cfl_s', 0.25) / (rel / dsz).max(), cfg.get('dtmax', 0.04))
        if t >= tout - 1e-9:
            KE, PE, vol = sol.energy(zeta, phi)
            if E0 is None:
                E0 = KE + PE
            gp = min_gap(z)
            frames.append(dict(t=round(t, 5), x=z.real.copy(), y=z.imag.copy(), a=a.copy(),
                               u=uu.real.copy(), v=uu.imag.copy(), phi=phi.copy()))
            if verbose:
                print(f"t={t:6.3f} N={N} dt={dt:.4f} E={(KE+PE)/E0:.5f} vol={vol:8.3f} ymax={z.imag.max():.3f} "
                      f"ymin={z.imag.min():.3f} gap={gp:.3f} dsmin={dsz.min():.3f} wall={time.time()-t0:.0f}s", flush=True)
            dzz = np.diff(z)
            steepness = np.abs(dzz.imag / np.abs(dzz)).max()
            tout += dto if steepness < 0.4 else cfg.get('dto_fine', 0.02)
            if gp < gap_stop:
                print("touchdown imminent, stop")
                break
        dt = min(dt, max(tout - t, 1e-4))
        zeta, phi = sol.step(zeta, phi, dt)
        t += dt
        step += 1
        st = cfg.get('filt', 36.0)
        e = np.arange(N) * L / N
        zeta = lowpass(zeta - e, st) + e
        phi = lowpass(phi, st).real
        if not np.isfinite(zeta).all():
            print("NaN, stop")
            break
    bot = reef.bottom()
    np.savez_compressed(out, L=L, h=reef.h, bx=bot.real, by=bot.imag,
                        t=np.array([f['t'] for f in frames]),
                        n=np.array([len(f['x']) for f in frames]),
                        x=np.concatenate([f['x'] for f in frames]),
                        y=np.concatenate([f['y'] for f in frames]),
                        a=np.concatenate([f['a'] for f in frames]),
                        u=np.concatenate([f['u'] for f in frames]),
                        v=np.concatenate([f['v'] for f in frames]))
    return frames, reef


def ap_fix(ap, zeta, L):
    """label was interpolated as (a - L e/2pi) vs old param; that is ambiguous after
    reparametrisation, so we carry d = a - Re(zeta) instead (see run)."""
    return ap + zeta.real
