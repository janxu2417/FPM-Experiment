import json
import math
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat


ROOT = Path(__file__).resolve().parents[2]
CODE_DIR = ROOT / "Code"
OUT_DIR = Path(__file__).resolve().parent


def mround(x):
    return np.floor(np.asarray(x) + 0.5).astype(int)


def center_matrix(a, m, n):
    out = np.zeros((m, n), dtype=a.dtype)
    r0 = (m - a.shape[0]) // 2
    c0 = (n - a.shape[1]) // 2
    out[r0:r0 + a.shape[0], c0:c0 + a.shape[1]] = a
    return out


def gaussian_blur_separable(img, sigma_px):
    if sigma_px <= 0:
        return img
    rad = max(1, int(math.ceil(sigma_px * 3)))
    ax = np.arange(-rad, rad + 1)
    kernel = np.exp(-0.5 * (ax / sigma_px) ** 2)
    kernel /= kernel.sum()
    out = np.apply_along_axis(lambda v: np.convolve(v, kernel, mode="same"), 0, img)
    out = np.apply_along_axis(lambda v: np.convolve(v, kernel, mode="same"), 1, out)
    return out


def measure_resolution_groups(img, spsize):
    freqs = [330, 300, 270, 240, 220]
    current_y_center = 100e-6
    m, n = img.shape
    x = (np.arange(n) - n / 2 + 0.5) * spsize
    y = (np.arange(m) - m / 2 + 0.5) * spsize
    rows = []

    for ii, f in enumerate(freqs):
        T = 1.0 / (f * 1000.0)
        x_start_v = -15e-6
        x_end_v = x_start_v + 4.5 * T
        y_start_v = current_y_center - 3.5 * T
        y_end_v = current_y_center + 3.5 * T

        x_start_h = x_end_v + T
        x_end_h = x_start_h + 7 * T
        y_start_h = y_start_v
        y_end_h = y_start_h + 4.5 * T

        c_v = line_contrast(img, x, y, (x_start_v, x_end_v, y_start_v, y_end_v), "v", T)
        c_h = line_contrast(img, x, y, (x_start_h, x_end_h, y_start_h, y_end_h), "h", T)

        rows.append(
            {
                "freq_lpmm": f,
                "contrast_v": float(c_v),
                "contrast_h": float(c_h),
                "contrast_avg": float(0.5 * (c_v + c_h)),
                "contrast_max": float(max(c_v, c_h)),
            }
        )
        if ii < len(freqs) - 1:
            T_next = 1.0 / (freqs[ii + 1] * 1000.0)
            current_y_center = current_y_center - 4.5 * T - 5.5 * T_next
    return rows


def line_contrast(img, x, y, rect, orient, T):
    x0, x1, y0, y1 = rect
    xmask = (x >= x0) & (x <= x1)
    ymask = (y >= y0) & (y <= y1)
    roi = img[np.ix_(ymask, xmask)]

    if orient == "v":
        profile = roi.mean(axis=0)
        coord = x[xmask]
        origin = x0
    else:
        profile = roi.mean(axis=1)
        coord = y[ymask]
        origin = y0

    phase = np.mod(coord - origin, T)
    bright = profile[phase < T / 2]
    dark = profile[phase >= T / 2]
    ib = bright.mean()
    idk = dark.mean()
    return abs((ib - idk) / max(ib + idk, np.finfo(float).eps))


def write_ppm(path, img, title):
    arr = np.clip(img, 0, 1)
    arr = np.kron(arr, np.ones((3, 3)))
    arr[:18, :] = 1.0
    x = 4
    for ch in title.encode("ascii", "ignore")[:24]:
        bits = [(ch >> b) & 1 for b in range(7, -1, -1)]
        for bit in bits:
            if x < arr.shape[1] - 2 and bit:
                arr[3:15, x:x + 2] = 0.0
            x += 2
        x += 2
    rgb = np.stack([arr, arr, arr], axis=-1)
    rgb = np.clip(rgb * 255, 0, 255).astype(np.uint8)
    with open(path, "wb") as f:
        f.write(f"P6\n{rgb.shape[1]} {rgb.shape[0]}\n255\n".encode("ascii"))
        f.write(rgb.tobytes())


def main():
    spsize = 1.38e-6
    pratio = 10
    psize = spsize / pratio
    m1 = n1 = 300
    m = m1 * pratio
    n = n1 * pratio
    wlength = 628e-9
    na_obj = 0.15
    z = 2e-6
    k0 = 2 * np.pi / wlength

    cfg = {
        "partial_coherence": {
            "enable_image_blur": True,
            "image_sigma_px": 0.8,
        },
        "pupil_mtf": {
            "enable": True,
            "rolloff_sigma": 0.85,
            "rolloff_power": 2.0,
        },
    }

    o_ideal = loadmat(CODE_DIR / "Raw_input" / "FP_ideal_data.mat")["O_ideal"]
    i_hr = center_matrix(o_ideal, m, n)
    i_hr = np.maximum(i_hr, 0.01)
    i_hr = i_hr / i_hr.max()
    cy, cx = i_hr.shape[0] // 2, i_hr.shape[1] // 2
    i_crop = i_hr[cy - m // 2:cy + m // 2, cx - n // 2:cx + n // 2]
    amp = np.sqrt(i_crop)
    phase = (1 - amp) * np.pi
    obj = amp * np.exp(1j * phase)
    obj_ft = np.fft.fftshift(np.fft.fft2(obj))

    kx, ky = np.meshgrid(
        (np.arange(n1) - n1 / 2) * (2 * np.pi / (n1 * spsize)),
        (np.arange(m1) - m1 / 2) * (2 * np.pi / (m1 * spsize)),
    )
    kval = np.sqrt(kx ** 2 + ky ** 2)
    kc = na_obj * k0
    rho = kval / max(kc, np.finfo(float).eps)
    pupil = (kval <= kc).astype(float)
    if cfg["pupil_mtf"]["enable"]:
        pupil *= np.exp(-((rho / cfg["pupil_mtf"]["rolloff_sigma"]) ** cfg["pupil_mtf"]["rolloff_power"]))
    pupil = pupil * np.exp(1j * kval ** 2 * z / (2 * k0))

    kxc = mround((n + 1) / 2)
    kyc = mround((m + 1) / 2)
    kyl = mround(kyc - (m1 - 1) / 2) - 1
    kyh = mround(kyc + (m1 - 1) / 2)
    kxl = mround(kxc - (n1 - 1) / 2) - 1
    kxh = mround(kxc + (n1 - 1) / 2)

    field = np.fft.ifft2(np.fft.ifftshift(obj_ft[kyl:kyh, kxl:kxh] * pupil))
    sim_center = np.abs(field)
    if cfg["partial_coherence"]["enable_image_blur"]:
        sim_center = gaussian_blur_separable(sim_center, cfg["partial_coherence"]["image_sigma_px"])
    sim_center = sim_center / sim_center.max()

    exp_center = loadmat(CODE_DIR / "FPM_RawData_3HDR.mat")["imlow_HDR"][:, :, 0].astype(float)
    exp_center = exp_center / exp_center.max()

    metrics_sim = measure_resolution_groups(sim_center, spsize)
    metrics_exp = measure_resolution_groups(exp_center, spsize)

    savemat(
        OUT_DIR / "nonideal_center_brightfield_result.mat",
        {
            "sim_center": sim_center,
            "exp_center": exp_center,
            "metrics_sim_freq_lpmm": np.array([r["freq_lpmm"] for r in metrics_sim]),
            "metrics_sim_contrast_avg": np.array([r["contrast_avg"] for r in metrics_sim]),
            "metrics_sim_contrast_max": np.array([r["contrast_max"] for r in metrics_sim]),
            "metrics_exp_contrast_avg": np.array([r["contrast_avg"] for r in metrics_exp]),
            "metrics_exp_contrast_max": np.array([r["contrast_max"] for r in metrics_exp]),
        },
    )

    with open(OUT_DIR / "nonideal_center_brightfield_metrics.json", "w", encoding="utf-8") as f:
        json.dump(
            {
                "config": cfg,
                "simulation": metrics_sim,
                "experiment": metrics_exp,
                "threshold": 0.03,
            },
            f,
            ensure_ascii=False,
            indent=2,
        )

    write_ppm(OUT_DIR / "nonideal_center_brightfield_sim.ppm", sim_center, "sim_nonideal")
    write_ppm(OUT_DIR / "nonideal_center_brightfield_exp.ppm", exp_center, "exp_center")


if __name__ == "__main__":
    main()
