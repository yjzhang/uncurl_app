import json
from multiprocessing.dummy import Process
import os
import uuid

from flask import Markup, render_template, request, redirect, send_from_directory, url_for, Blueprint, current_app
from werkzeug.utils import secure_filename

import numpy as np
import uncurl

from .cache import cache

from .generate_analysis import generate_uncurl_analysis, get_progress
from .data_stats import Summary

views = Blueprint('views', __name__, template_folder='templates')

# TODO: allow this to upload multiple files with multiple genes
def load_upload_data(request_files, request_form, path=None):
    # TODO: should we merge the datasets here... or merge them in a downstream step?
    data_paths = []
    gene_paths = []
    shapes = []
    output_filenames = []
    data_names = []
    for key, f in request_files.items():
        if 'fileinput' in key:
            input_id = key.split('-')[1]
            output_filename = secure_filename(f.filename)
            if output_filename == '':
                return
            output_filenames.append(output_filename)
            # allow for mtx input data
            input_type = request_form['inputtype-{0}'.format(input_id)]
            if output_filename.endswith('.mtx.gz') or output_filename.endswith('.mtx'):
                input_type = 'sparse'
            if input_type == 'dense':
                data_filename = 'data_{0}.txt'.format(input_id)
                if output_filename.endswith('.gz'):
                    data_filename = 'data_{0}.txt.gz'.format(input_id)
                data_path = os.path.join(path, data_filename)
                f.save(data_path)
            elif input_type == 'sparse':
                data_filename = 'data_{0}.mtx'.format(input_id)
                if output_filename.endswith('.mtx.gz'):
                    data_filename = 'data_{0}.mtx.gz'.format(input_id)
                data_path = os.path.join(path, data_filename)
                f.save(data_path)
            shape = request_form['data_shape-{0}'.format(input_id)]
            shapes.append(shape)
            # TODO: try to find gene names
            data_paths.append(data_path)
            try:
                gene_file = request_files['genenames-{0}'.format(input_id)]
                print(gene_file)
                gene_output_file = load_gene_names(gene_file, path, input_id)
                gene_paths.append(gene_output_file)
            except:
                gene_paths.append(None)
            data_name = request_form['data_name-{0}'.format(input_id)]
            if len(data_name) == 0:
                data_name = str(input_id)
            data_names.append(data_name)
    init = None
    return data_paths, gene_paths, data_names, init, shapes


def load_gene_names(f, path=None, index=1):
    gene_filename = secure_filename(f.filename)
    gene_output_path = os.path.join(path, 'gene_names_{0}.txt'.format(index))
    if gene_filename.endswith('genes.csv'):
        import pandas as pd
        data = pd.read_csv(f)
        try:
            gene_names = data['gene_name']
            gene_names.to_csv(gene_output_path, header=None, index=None)
        except Exception as e:
            print(e)
            f.seek(0)
            f.save(gene_output_path)
    elif gene_filename.endswith('features.tsv'):
        import pandas as pd
        data = pd.read_csv(f, sep='\t', header=None)
        try:
            gene_names = data[1]
            gene_names.to_csv(gene_output_path, header=None, index=None)
        except Exception as e:
            print(e)
            f.seek(0)
            f.save(gene_output_path)
    else:
        f.save(gene_output_path)
    return gene_output_path

@views.route('/help')
def help():
    return render_template('help.html')

@views.route('/state_estimation')
@cache.cached()
def state_estimation():
    return render_template('state_estimation.html')

@views.route('/state_estimation/input', methods=['POST'])
def state_estimation_input():
    user_id = str(uuid.uuid4())
    if 'username' in request.form:
        if len(request.form['username']) > 0:
            # make username a safe string
            keep_chars = set(['-', '_', ' '])
            username = request.form['username'].strip()[:25]
            username = ''.join([c for c in username if c.isalnum() or (c in keep_chars)])
            user_id = user_id + '-' + username
    base_path = os.path.join(current_app.config['USER_DATA_DIR'], user_id)
    os.makedirs(base_path)
    # save request.form
    with open(os.path.join(base_path, 'inputs.json'), 'w') as f:
        f.write(json.dumps(request.form))
    # TODO: if file is large, start a new thread. otherwise just
    # run the thing
    request_file = request.files
    request_form = request.form
    data_paths, gene_paths, output_filenames, init, shapes = load_upload_data(request_file, request_form, base_path)
    # TODO: deal with init
    P = Process(target=state_estimation_preproc, args=(user_id, base_path, data_paths, gene_paths, output_filenames, init,
        shapes))
    P.start()
    #state_estimation_preproc(user_id, path)
    return redirect(url_for('views.state_estimation_result', user_id=user_id))

@views.route('/state_estimation/results/<user_id>/start', methods=['POST'])
def state_estimation_start(user_id):
    """
    Actually start the process of state estimation.

    This saves a file called 'params.json' in /tmp/uncurl/<user_id>
    containing all parameters used in state estimation.
    """
    path = os.path.join(current_app.config['USER_DATA_DIR'], user_id)
    gene_names_file = os.path.join(path, 'gene_names.txt')
    if not os.path.exists(gene_names_file):
        gene_names_file = None
    # TODO: deal with init here - make note if it's qualitative or
    # quantitative
    # run qualNorm???
    init_path= os.path.join(path, 'init.txt')
    if not os.path.exists(init_path):
        init_path = None
    # load json params
    with open(os.path.join(path, 'preprocess.json')) as f:
        preprocess = json.load(f)
    for key in request.form.keys():
        preprocess[key] = request.form[key]
    # params.json contains all input parameters to the state estimation, as well as all stats from preprocess.json.
    with open(os.path.join(path, 'params.json'), 'w') as f:
        json.dump(preprocess, f)
    P = Process(target=state_estimation_thread, args=(user_id, gene_names_file, init_path, path, preprocess, current_app.config.copy()))
    P.start()
    return redirect(url_for('views.state_estimation_result', user_id=user_id))


@views.route('/state_estimation/results/<user_id>/')
def state_estimation_result(user_id):
    path = os.path.join(current_app.config['USER_DATA_DIR'], user_id)
    if os.path.exists(os.path.join(path, 'sc_analysis.json')):
        return redirect(url_for('interaction_views.view_plots', user_id=user_id))
    elif os.path.exists(os.path.join(path, 'preprocess.json')):
        uncurl_is_running = os.path.exists(os.path.join(path, 'submitted'))
        current_task = 'None'
        time_remaining = 'Unknown'
        with open(os.path.join(path, 'preprocess.json')) as f:
            preprocess = json.load(f)
        try:
            with open(os.path.join(path, 'read_count_hist_data.json')) as f:
                read_count_hist_data = f.read()
            with open(os.path.join(path, 'gene_count_hist_data.json')) as f:
                gene_count_hist_data = f.read()
            with open(os.path.join(path, 'gene_mean_hist_data.json')) as f:
                gene_mean_hist_data = f.read()
        except:
            summary = Summary(None, None, base_path=path)
            read_count_hist_data, gene_count_hist_data, gene_mean_hist_data = summary.generate_plotly_jsons()
        if uncurl_is_running:
            # get running time information (highly approximate)
            current_task, time_remaining = get_progress(path)
            # update with actual input parameters
            with open(os.path.join(path, 'params.json')) as f:
                preprocess.update(json.load(f))
        uncurl_has_error = False
        if time_remaining == 'error':
            # TODO: if there is an error here...
            uncurl_has_error = True
        return render_template('state_estimation_user.html',
                user_id=user_id, has_preview=True,
                uncurl_is_done=False,
                uncurl_is_running=uncurl_is_running,
                uncurl_has_error=uncurl_has_error,
                read_count_hist_data=read_count_hist_data,
                gene_count_hist_data=gene_count_hist_data,
                gene_mean_hist_data=gene_mean_hist_data,
                current_task=current_task,
                time_remaining=time_remaining,
                **preprocess)
    elif os.path.exists(os.path.join(path, 'error.txt')):
        error_txt = ''
        with open(os.path.join(path, 'error.txt')) as f:
            error_txt = f.read()
        return error(error_txt, 404)
    else:
        return render_template('state_estimation_user.html',
                user_id=user_id, uncurl_is_running=False,
                uncurl_is_done=False,
                has_result=False)

# this gzips the directory and returns a download
@views.route('/<x>/results/<user_id>/download_all')
def state_estimation_download_all(x, user_id):
    if x!='test':
        path = os.path.join(current_app.config['USER_DATA_DIR'], user_id)
        if not os.path.exists(path):
            path = os.path.join(current_app.config['SECONDARY_USER_DATA_DIR'], user_id)
            if not os.path.exists(path):
                return error('Data not found', 404)
    else:
        path = os.path.join(current_app.config['TEST_DATA_DIR'], user_id)
    filename = user_id + '.tar.gz'
    output_filename = os.path.join(current_app.config['USER_DATA_DIR'], filename)
    create_tar = True
    # update tarball if path is newer than output_filename
    if os.path.exists(output_filename):
        tar_mtime = os.stat(output_filename.st_mtime)
        create_tar = False
        for base, dirs, files in os.walk(path):
            for f in files:
                if os.stat(os.path.join(base, f)).st_mtime > tar_mtime:
                    create_tar = True
                    break
    if create_tar:
        import subprocess
        subprocess.call(['tar', '-czf', output_filename, path])
    print('download_all', path, filename)
    return send_from_directory(current_app.config['USER_DATA_DIR'], filename)

@views.route('/<x>/results/<user_id>/<filename>')
def state_estimation_file(x, user_id, filename):
    if x != 'test':
        path = os.path.join(current_app.config['USER_DATA_DIR'], user_id)
    else:
        path = os.path.join(current_app.config['TEST_DATA_DIR'], user_id)
    print('download: ', path)
    return send_from_directory(path, filename)

@views.route('/<x>/results/<user_id>/data_download')
def data_download(x, user_id):
    if x!='test':
        path = os.path.join(current_app.config['USER_DATA_DIR'], user_id)
    else:
        path = os.path.join(current_app.config['TEST_DATA_DIR'], user_id)
    files = os.listdir(path)
    files.sort()
    return render_template('data_download.html',
            user_id=user_id,
            test_or_user=x,
            files=files)

def state_estimation_preproc(user_id, base_path, data_paths, gene_paths, output_filenames, init=None,
        shapes=['gene_cell']):
    # TODO: update for multiple data/genes
    """
    Preprocessing for state estimation - generates summary statistics,
    etc...
    """
    # TODO: combine datasets
    if base_path is None:
        base_path = os.path.join(current_app.config['USER_DATA_DIR'], user_id)
    try:
        summary = Summary(data_paths, gene_paths, base_path, shapes=shapes, dataset_names=output_filenames)
        read_count_hist_data, gene_count_hist_data, gene_mean_hist_data = summary.load_plotly_json()
        summary.preprocessing_params()
    except:
        import traceback
        text = traceback.format_exc()
        with open(os.path.join(base_path, 'error.txt'), 'w') as f:
            f.write(text)

def state_estimation_preproc_simple(user_id, base_path, data_path):
    """
    Preprocessing, assuming that the data has already been merged.
    """
    # TODO: combine datasets
    if base_path is None:
        base_path = os.path.join(current_app.config['USER_DATA_DIR'], user_id)
    try:
        summary = Summary(None, None, base_path)
        read_count_hist_data, gene_count_hist_data, gene_mean_hist_data = summary.load_plotly_json()
        summary.preprocessing_params()
    except:
        import traceback
        text = traceback.format_exc()
        with open(os.path.join(base_path, 'error.txt'), 'w') as f:
            f.write(text)



def state_estimation_thread(user_id, gene_names=None, init_path=None, path=None, preprocess=None, config=None):
    """
    Uses a new process to do state estimation. Assumes that the input data is already saved in a directory named /tmp/uncurl/<user_id>/.

    Args:
        user_id (str)
        gene_names (str or array, optional): path to a list of gene names, or an array of gene names. Default: None
        init_path (str, optional): path to txt matrix of shape (genes, k). Default: None.
        path (str, optional): Path where data and results are saved.
        preprocess (dict): dict containing additional parameters: min_reads, max_reads, normalize, is_sparse, is_gz, disttype, genes_frac, cell_frac, vismethod, baseline_vismethod
        config (dict): current_app.config
    """
    if path is None:
        path = os.path.join(config['USER_DATA_DIR'], user_id)
    # get correct data names
    data = os.path.join(path, 'data.mtx')
    if preprocess['is_gz']:
        data += '.gz'
    if not os.path.exists(data):
        data = os.path.join(path, 'data.txt')
        if preprocess['is_gz']:
            data += '.gz'
    # TODO: it's really confusing where the param names come from in
    # the preprocess dict - they come from the input names in
    # state_estimation_user.html.
    dist_type = preprocess['disttype']
    if dist_type=='Poisson':
        pass
    elif dist_type=='Negative binomial':
        dist_type = 'nb'
    elif dist_type == 'Log-Normal':
        dist_type = 'lognorm'
    uncurl_args = config['UNCURL_ARGS']
    if dist_type != 'Poisson':
        uncurl_args = config['NMF_ARGS']
    if dist_type == 'Poisson':
        uncurl_args['write_progress_file'] = os.path.join(path, 'progress.txt')
    uncurl_args['dist'] = dist_type
    # TODO: deal with init
    if init_path is not None:
        pass
    # params.json is saved in path, so it does not need to be passed.
    generate_uncurl_analysis(data, path,
            **uncurl_args)


@views.route('/qual2quant')
def qual2quant():
    return render_template('qual2quant.html')

@views.route('/qual2quant/input')
def qual2quant_input():
    if 'fileinput' not in request.files or 'qualinput' not in request.files:
        return error('Missing data input', 400)
    cell_file = request.files['fileinput']
    qual_file = request.files['qualinput']
    cell_data = np.loadtxt(cell_file)
    qual_data = np.loadtxt(qual_file)
    user_id = str(uuid.uuid4())
    P = Process(target=qual2quant_thread, args=(cell_data, qual_data, user_id))
    P.start()
    return redirect(url_for('qual2quant_result', user_id=user_id))

def qual2quant_thread(data, qual, user_id):
    centers = uncurl.qualNorm(data, qual)
    path = os.path.join(current_app.config['USER_DATA_DIR'], user_id)
    with open(os.join(path, 'qual2quant_centers.txt'), 'w') as f:
        np.savetxt(f, centers)

@views.route('/qual2quant/results/<user_id>')
def qual2quant_result(user_id):
    if os.path.exists(os.path.join(current_app.config['USER_DATA_DIR'], user_id, 'qual2quant_centers.txt')):
       return render_template('qual2quant_user.html',
                user_id=user_id, has_result=True,
                visualization=None)
    else:
        return render_template('qual2quant_user.html',
                user_id=user_id, has_result=False)


def error(msg, code):
    return render_template('error.html', msg=msg), code
