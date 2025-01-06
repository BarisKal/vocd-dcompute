import random
import math
from typing import List
import matplotlib.pyplot as plt


def find_min_d(d_av, nt_dict):
    """
    Find the value of d that minimizes the least squares difference
    between the t-equation and a set of <n,t> values.

    Parameters:
    d_av (float): The seeded d value.
    nt_dict (list of dict): List of dictionaries containing 'N' and 'T' values.

    Returns:
    tuple: The value of D for which the least squares difference is a minimum,
           and the minimum least squares difference.
    """
    index = 0

    stepsize = 0.001
    diff = d_least_sq(d_av, nt_dict) - d_least_sq(d_av - stepsize, nt_dict)

    if diff > 0:
        k = -1
    elif diff < 0:
        k = 1
    else:
        # if k == 0
        return d_av, d_least_sq(d_av, nt_dict)

    prev_d_least_sq = d_av
    d = d_av
    while d > 0 and d < 2 * d_av:
        print(index, d)
        next_ls = d_least_sq(d, nt_dict)
        if prev_d_least_sq < next_ls:
            print("Breaking at", index, d, d_av, prev_d_least_sq, next_ls)
            break

        prev_d_least_sq = next_ls
        index += 1
        d += k * 0.001

    print(d_av, d)

    return d, prev_d_least_sq


def d_least_sq(d, nt_dict):
    """
    Calculate the least squares difference for a given d value.

    Parameters:
    d (float): The d value.
    nt_dict (list): List of dictionaries containing 'N' and 'T' values.

    Returns:
    float: The least squares difference.
    """
    return sum((nt["TTR"] - ttr_eqn(d, nt["N"])) ** 2 for nt in nt_dict)


def ttr_eqn(d, n):
    """
    Calculate the ttr value for a given d and n.

    Parameters:
    d (float): The d value.
    n (int): The n value.

    Returns:
    float: The ttr value.
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
            TTR, SD = average_ttr(token_seq, tkns_in_seq, N, no_samples)
            nt_dict.append({"N": N, "TTR": TTR, "SD": SD, "D": d_eqn(N, TTR)})

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

    # Sample text from \clan\examples\transcripts\ne32\68.cha
    test_string = """what be in these box oh what be in it yeah ahhah that be just like mine at home yeah let us read it like mine yeah Frannie haha yeah yeah uh yeah I there ahhah yeah here be one for you I have Cookie_Monster my name be Cookie_Monster and I eat cookie nope I do not yeah except they go up and down hm yup moo yeah yeah oh oh we can draw thing on paper oh there we get orange orange no a mouth that could be the lady have the on how about you draw the lady with the hat no make that big big that big here be some blue that big oh yeah I want to draw her lip she will sing with no you make the lip right here no right here yeah she be sing with open yeah yeah no they be all do yeah I want to make that lady she be sing on tv I want to do the lady too that one no you do it no you do it no I want to do that one this be so heavy oh that look like mine too nope in here yup it do okay who be home in here Maude right here okay she want to go in the chair a chair on the table table table table table table that be my seat you guy I want to sit in them nope a car that be mine no he where be his seat he be get angry too woofwoof where be my chair I want vroom there can not reach the chair there who take my chair in that chair you do not believe it I do not believe it where be my seat he can take his own chair if he want the green chair let me have the green chair woofwoof we think dog do not like person to go in the chair hm knock knock hm nope maybe I could go back in that box he want to go in that box yeah Harold hide you get those you get Maude please please yeah no he be hide he be angry no I do not want to come out I tire I want to sleep yeah no I have my dinner this kid could come to the table"""

    splitted_text = test_string.split(" ")
    d_optimum_avg = d_compute(splitted_text, len(splitted_text), 35, 50, 1, 100)
    print("D optimum average: " + str(d_optimum_avg))
