import random
import math
from typing import List
from dataclasses import dataclass

test_string = """
Once there was a kind man. He had a wife and a daughter. His wife loved him very much, and so did his daughter. One day the man's wife died. He was very, very sad.
The man wanted a new wife. His next wife was not kind at all. She was very cruel. She had two ugly daughters who were also very mean.
Once there was a kind man. He had a wife and a daughter. His wife loved him very much, and so did his daughter. One day the man's wife died. He was very, very sad.
The man wanted a new wife. His next wife was not kind at all. She was very cruel. She had two ugly daughters who were also very mean.
Once there was a kind man. He had a very loving wife. He had a young daughter who was also just like her mother. One day the man's wife died, and he was very sad.
The wife did not like her new daughter, who was kind and pretty. So she made her work hard. The girl scrubbed dishes. She scrubbed floors. She cleaned the fireplaces. Her sisters made fun of her and called her 'ash girl', or Cinderella.
Cinderella's sisters had fine rooms with soft beds. But Cinderella had a cold room in the attic, and her bed was made of straw.
"""


def find_min_d(d_av, nt_tuples):
    """
    Find the value of d that minimizes the least squares difference
    between the t-equation and a set of <n,t> values.

    Parameters:
    d_av (float): The seeded d value.
    nt_tuples (list): List of dictionaries containing 'N' and 'T' values.

    Returns:
    tuple: The value of D for which the least squares difference is a minimum,
           and the minimum least squares difference.
    """
    stepsize = 0.001
    diff = d_least_sq(d_av, nt_tuples) - d_least_sq(d_av - stepsize, nt_tuples)

    if diff > 0:
        k = -1
    elif diff < 0:
        k = 1
    else:
        # if k == 0
        return d_av, d_least_sq(d_av, nt_tuples)

    prev_d_least_sq = d_av
    for d in [d_av + k * stepsize * i for i in range(int(2 * d_av / stepsize))]:
        next_ls = d_least_sq(d, nt_tuples)
        if prev_d_least_sq < next_ls:
            break
        prev_d_least_sq = next_ls

    return d, prev_d_least_sq


def d_least_sq(d, nt_tuples):
    """
    Calculate the least squares difference for a given d value.

    Parameters:
    d (float): The d value.
    nt_tuples (list): List of dictionaries containing 'N' and 'T' values.

    Returns:
    float: The least squares difference.
    """
    return sum((nt["T"] - t_eqn(d, nt["N"])) ** 2 for nt in nt_tuples)


def t_eqn(d, n):
    """
    Calculate the t value for a given d and n.

    Parameters:
    d (float): The d value.
    n (int): The n value.

    Returns:
    float: The t value.
    """
    x = math.sqrt(1 + (2 * n) / d) - 1
    return (d / n) * x


def d_eqn(n: int, ttr: float) -> float:
    """
    Calculate the D value based on the given parameters.
    This function computes the D value using the formula:
    D = 0.5 * ((n * ttr^2) / (1 - ttr))
    Parameters:
    n (int): The number of tokens.
    ttr (float): The type-token ratio.
    Returns:
    float: The computed D value. If the type-token ratio is 1, the function returns 0 to avoid division by zero.
    """
    d = 0
    if (1 - ttr) == 0:
        return 0

    d = 0.5 * ((n * ttr * ttr) / (1 - ttr))
    return d


def average_ttr(token_seq, tkns_in_seq, segment_size, no_samples):
    """
    Computes the average TTR (type-token ratio) and standard deviation for the given segment size.

    Parameters:
    token_seq (list): List of tokens.
    tkns_in_seq (int): Total number of tokens in the sequence.
    segment_size (int): Size of each segment.
    no_samples (int): Number of samples.
    flags (int): Flags indicating the mode of operation.

    Returns:
    tuple: Average TTR and standard deviation.
    """
    no_segments = 100  # Taken from D_samples
    ttr_list = []

    for _ in range(no_segments):
        # Sampling without replacement; segment_size will be 35 to 50
        segment = random.sample(token_seq, segment_size)

        types = len(set(segment))
        ttr_list.append(types / segment_size)

    avg_ttr = sum(ttr_list) / len(ttr_list)
    std_dev = math.sqrt(sum((x - avg_ttr) ** 2 for x in ttr_list) / len(ttr_list))
    return avg_ttr, std_dev


def d_compute(
    token_seq: List[str],
    tkns_in_seq: int,
    from_d: int,
    to_d: int,
    incr: int,
    no_samples: int,
) -> int:
    """
    This function computes the values of d based on the provided sequence of tokens and other parameters.
    If the flag is RANDOM, it computes three values using a defined set of samples selected at random,
    of a particular size. If the flag is not RANDOM, it computes a single value of d.
    Parameters:
    token_seq (List[str]): The sequence of tokens.
    tkns_in_seq (int): The number of tokens in the sequence.
    from_d (int): The lower value of d.
    to_d (int): The upper value of d.
    incr (int): The step size of d.
    no_samples (int): The number of segments to be used in the calculation of the average value of TTR.

    Returns:
    int: The computed average minimum value of d.
    """
    no_NTvalues = 0
    no_trials = 3
    dmin_trials = []

    for j in range(no_trials):

        # calculate the number of <n,t> values to be returned
        no_NTvalues = ((to_d - from_d) / incr) + 1
        nt_dict = []

        for N in range(from_d, to_d + 1, incr):
            T, SD = average_ttr(token_seq, tkns_in_seq, N, no_samples)
            nt_dict.append({"N": N, "T": T, "SD": SD, "D": d_eqn(N, T)})

        # Count in nt_dict if D is 0 and discard it while calculating the average of D
        # Usually D shouldn't be exactly 0 but in the original implementation it was defined. Probably to catch edge cases.
        discard_d = sum(1 for nt in nt_dict if nt["D"] == 0)

        d_av = sum(nt["D"] for nt in nt_dict if nt["D"] != 0) / (
            len(nt_dict) - discard_d
        )

        d_std = sum(math.pow((nt["D"] - d_av), 2) for nt in nt_dict)
        d_std = math.sqrt(
            d_std / (no_NTvalues - discard_d - 1)
        )  # Bessel's correction not available in the original implementation

        d_min, min_ls_val = find_min_d(d_av, nt_dict)
        dmin_trials.append(d_min)

    dmin_av = sum(dmin_trials) / len(dmin_trials)

    return dmin_av


if __name__ == "__main__":
    mydmin = d_compute(test_string.split(" "), len(test_string.split()), 35, 50, 1, 100)
    print(mydmin)
