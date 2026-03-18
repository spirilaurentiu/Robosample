import re
import numpy as np

def parse_prmtop_numpy(prmtop_file):
    FORMAT_RE_PATTERN = re.compile(r"(\d+)\(?([a-zA-Z]+)(\d+)\.?(\d*)\)?")

    flags = []
    raw_data = {}
    raw_format = {}
    prmtop_version = None

    with open(prmtop_file, 'r') as f:
        lines = [line.rstrip('\n') for line in f]

    for line in lines:
        if not line:
            continue
        if line.startswith('%'):
            if line.startswith('%VERSION'):
                _, prmtop_version = line.split(None, 1)
            elif line.startswith('%FLAG'):
                _, flag = line.split(None, 1)
                flag = flag.strip()
                flags.append(flag)
                raw_data[flag] = [] 
            elif line.startswith('%FORMAT'):
                fmt_line = line[line.index('(')+1 : line.index(')')]
                m = FORMAT_RE_PATTERN.search(fmt_line)
                if m:
                    raw_format[flags[-1]] = (fmt_line, int(m.group(1)), m.group(2),int(m.group(3)), m.group(4))
                else:
                    raw_format[flags[-1]] = (fmt_line, 1, 'a', 80, '')
            continue

        # Non-comment, non-flag lines -> data
        flag = flags[-1]
        fmt, num_items, item_type, i_length, item_prec = raw_format[flag]

        if flag == 'TITLE' and not raw_data[flag]:
            raw_data[flag] = [line]
            continue

        # Vectorized chunking
        arr = np.frombuffer(line.encode('utf-8'), dtype='S1')
        n_chunks = len(arr) // i_length + (len(arr) % i_length > 0)
        chunks = [arr[i*i_length:(i+1)*i_length].tobytes().decode('utf-8') for i in range(n_chunks)]
        items = [c.strip() for c in chunks if c.strip()]

        # Convert according to type
        if item_type.upper() == 'A':
            raw_data[flag].extend(items)
        elif item_type.upper() == 'I':
            raw_data[flag].extend(np.array(items, dtype=np.int64))
        elif item_type.upper() in ('E', 'F', 'D'):
            # Fortran-style floats, parse as float64
            raw_data[flag].extend(np.array([float(x.replace('D', 'E')) for x in items], dtype=np.float64))
        else:
            # fallback as string
            raw_data[flag].extend(items)

    # Convert numeric lists to np.array for consistency
    for flag in raw_data:
        if isinstance(raw_data[flag], list) and raw_data[flag]:
            first_item = raw_data[flag][0]
            if isinstance(first_item, (int, np.integer)):
                raw_data[flag] = np.array(raw_data[flag], dtype=np.int64)
            elif isinstance(first_item, (float, np.floating)):
                raw_data[flag] = np.array(raw_data[flag], dtype=np.float64)

    # Per AMBER prmtop convention, atomic charges are stored multiplied by 18.2223
    # We divide to recover physical charges in units of the proton/electron charge.
    raw_data['CHARGE'] /= 18.2223

    chamber_style = 'CTITLE' in flags
    return {
        'version': prmtop_version,
        'flags': flags,
        'raw_data': raw_data,
        'raw_format': raw_format,
        'chamber': chamber_style
    }

def has_nbfix_fast(nb_indices: np.ndarray, num_types: int, acoef: np.ndarray, bcoef: np.ndarray) -> bool:
    nb_indices = (
        np.array(nb_indices)
        .reshape(num_types, num_types) - 1
    )

    diag_idx = nb_indices.diagonal()
    A_ii = acoef[diag_idx]
    B_ii = bcoef[diag_idx]

    with np.errstate(divide='ignore', invalid='ignore'):
        rmin = (2 * A_ii / B_ii) ** (1/6)
        ei = 0.25 * B_ii**2 / A_ii

    ri = np.where(np.isfinite(rmin), rmin / 2.0, 0.0)
    ei = np.where(np.isfinite(ei), ei, 0.0)

    expected_R = ri[:, None] + ri[None, :]
    expected_E = np.sqrt(ei[:, None] * ei[None, :])

    mask = nb_indices >= 0

    actual_A = np.zeros((num_types, num_types))
    actual_B = np.zeros((num_types, num_types))
    actual_A[mask] = acoef[nb_indices[mask]]
    actual_B[mask] = bcoef[nb_indices[mask]]

    zero_mask = (actual_A == 0) | (actual_B == 0)
    bad_zero = zero_mask & (
        (actual_A != 0) |
        (actual_B != 0) |
        ((expected_E != 0) & (expected_R != 0))
    )

    if np.any(bad_zero & mask):
        return True

    calc_A = expected_E * expected_R**12
    calc_B = 2 * expected_E * expected_R**6

    bad_A = np.abs((actual_A - calc_A) / actual_A) > 1e-6
    bad_B = np.abs((actual_B - calc_B) / actual_B) > 1e-6

    return np.any((bad_A | bad_B) & mask)
