# Determine the Variance of the estimation (how much f(x) varies in the domain of x)
# this is an experiment: I had hypothesized the error would decrease if the same x value were sampled for the ave of squares
# and the sq of averages, however it resulted in a higher error.

def get_crude_mc_variance(num_samples=10000):
    """
    This function returns the variance for the Crude Monte Carlo
    """
    int_max = 5  # the max of our integration range

    # find the average of squares & sq of ave
    running_total_1 = 0
    running_total_2 = 0
    for i in range(num_samples):
        x = get_random_number(0, int_max)
        running_total_1 += f_of_x(x)**2
        running_total_2 += f_of_x(x)
    sum_of_sqs = (int_max * running_total_1 / num_samples)
    sq_ave = (int_max * running_total_2 / num_samples) ** 2

    return sum_of_sqs - sq_ave