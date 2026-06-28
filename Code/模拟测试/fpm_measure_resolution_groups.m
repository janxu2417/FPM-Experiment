function metrics = fpm_measure_resolution_groups(img, spsize)
% Measure contrast on the synthetic 330/300/270/240/220 lp/mm target.

freqs = [330, 300, 270, 240, 220];
current_y_center = 100e-6;

[m, n] = size(img);
x = ((1:n) - (n + 1) / 2) * spsize;
y = ((1:m) - (m + 1) / 2) * spsize;

metrics = struct();
for ii = 1:numel(freqs)
    f = freqs(ii);
    T = 1 / (f * 1000);

    x_start_v = -15e-6;
    x_end_v = x_start_v + 4.5 * T;
    y_start_v = current_y_center - 3.5 * T;
    y_end_v = current_y_center + 3.5 * T;

    x_start_h = x_end_v + T;
    x_end_h = x_start_h + 7 * T;
    y_start_h = y_start_v;
    y_end_h = y_start_h + 4.5 * T;

    c_v = local_line_contrast(img, x, y, [x_start_v, x_end_v, y_start_v, y_end_v], 'v', T);
    c_h = local_line_contrast(img, x, y, [x_start_h, x_end_h, y_start_h, y_end_h], 'h', T);

    metrics(ii).freq_lpmm = f;
    metrics(ii).contrast_v = c_v;
    metrics(ii).contrast_h = c_h;
    metrics(ii).contrast_avg = 0.5 * (c_v + c_h);
    metrics(ii).contrast_max = max(c_v, c_h);

    if ii < numel(freqs)
        T_next = 1 / (freqs(ii + 1) * 1000);
        current_y_center = current_y_center - 4.5 * T - 5.5 * T_next;
    end
end
end

function c = local_line_contrast(img, x, y, rect, orient, T)
x0 = rect(1); x1 = rect(2);
y0 = rect(3); y1 = rect(4);

xmask = (x >= x0) & (x <= x1);
ymask = (y >= y0) & (y <= y1);
roi = img(ymask, xmask);

if orient == 'v'
    profile = mean(roi, 1);
    coord = x(xmask);
    origin = x0;
else
    profile = mean(roi, 2).';
    coord = y(ymask);
    origin = y0;
end

phase = mod(coord - origin, T);
bright = profile(phase < T / 2);
dark = profile(phase >= T / 2);

Ib = mean(bright);
Id = mean(dark);
c = abs((Ib - Id) / max(Ib + Id, eps));
end
