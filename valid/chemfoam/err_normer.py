import numpy as np
import os

script_dir = os.path.dirname(os.path.abspath(__file__))


def path(file):
    return os.path.join(script_dir, file)


of_err_mean = None
of_err_inf = None
ai_err_mean = None
ai_err_inf = None
n_total = 0


def update_err(err, norm, n_samples, new):
    if norm == 'mean':
        # undo the square root and divide by the number of samples to get the mean
        new = new**2 / n_samples
    if err is None:
        err = new
    elif norm == 'mean':
        # add the squared error for later computation
        err += new
    else:
        assert norm == 'inf'
        # take the inf norm of the two
        err = np.maximum(np.abs(err), np.abs(new))
    return err


files = [path(x) for x in os.listdir(script_dir) if os.path.isfile(path(x))
         and x.endswith('.npz')]
for file in files:
    err = np.load(file)
    n_samp = err['n_samples']

    # and update
    of_err_mean = update_err(of_err_mean, 'mean', n_samp, err['of_err_mean'])
    of_err_inf = update_err(of_err_mean, 'inf', n_samp, err['of_err_inf'])

    ai_err_mean = update_err(ai_err_mean, 'mean', n_samp, err['ai_err_mean'])
    ai_err_inf = update_err(ai_err_inf, 'inf', n_samp, err['ai_err_inf'])

    n_total += n_samp

if len(files):
    # compute final mean error
    of_mean_comp = np.sqrt(of_err_mean / n_total)
    ai_mean_comp = np.sqrt(ai_err_mean / n_total)

    # and print component
    print('OpenFOAM mean (component) error:', np.array_str(of_mean_comp))
    print('AccelerInt mean (component) error:', np.array_str(ai_mean_comp))

    print('OpenFOAM max (component) error:', np.array_str(of_err_inf))
    print('AccelerInt max (component) error:', np.array_str(ai_err_inf))

    # and print general
    of_mean = np.linalg.norm(np.sqrt(of_err_mean / (n_total * len(files))), ord=2)
    ai_mean = np.linalg.norm(np.sqrt(ai_err_mean / (n_total * len(files))), ord=2)
    of_max = np.linalg.norm(of_err_inf, ord=np.inf, axis=0)
    ai_max = np.linalg.norm(ai_err_inf, ord=np.inf, axis=0)

    print('OpenFOAM mean error:', of_mean)
    print('AccelerInt mean error:', ai_mean)

    print('OpenFOAM max error:', of_max)
    print('AccelerInt max error:', ai_max)
