# This generates uncurl analyses as static HTML files.
# inputs: raw data (matrix or file), M (or file), W (or file), reduced_data (2d for vis)
import numpy as np
import scipy.io
from jinja2 import Environment, PackageLoader, select_autoescape

from bokeh.plotting import figure
from bokeh.embed import components
from bokeh.palettes import Accent8
from bokeh.models import CustomJS, ColumnDataSource, Div, Button


def display_event(div, attributes=[], style = 'float:left;clear:left;font_size=0.5pt'):
    """
    Build a suitable CustomJS to display the current event in the div model.
    Source: https://bokeh.pydata.org/en/latest/docs/user_guide/interaction/callbacks.html
    """
    # TODO: display given updates instead
    return CustomJS(args=dict(div=div), code="""
        var attrs = {0}; var args = [];
        for (var i=0; i<attrs.length; i++ ) {
        args.push(attrs[i] + '=' + Number(cb_obj[attrs[i]]).toFixed(2));
        }
        var line = "<span style={1}><b>" + cb_obj.event_name + "</b>(" + args.join(", ") + ")</span>\\n";
        var text = div.text.concat(line);
        var lines = text.split("\\n")
        if ( lines.length > 35 ) { lines.shift(); }
        div.text = lines.join("\\n");
        """.format(attributes, style))

def generate_analysis(data, m, w, reduced_data,
        gene_names=None,
        output_file='out.html',
        data_type='dense'):
    """
    Args:
        data (array): unprocessed data, shape (genes, cells).
        m (array): uncurl output, shape (genes, k).
        w (array): uncurl output, shape (k, cells).
        reduced_data (array): 2 x n array of data reduced either through tsne
            or something else
        gene_names (array): 1d array of strings, length=data.shape[0].
    """
    if isinstance(data, str):
        if data_type == 'dense':
            data = np.loadtxt(data)
        elif data_type == 'sparse':
            data = scipy.io.mmread(data)
    if isinstance(m, str):
        m = np.loadtxt(m)
    if isinstance(w, str):
        w = np.loadtxt(w)
    if isinstance(reduced_data, str):
        reduced_data = np.loadtxt(reduced_data)
    clusters = w.argmax(0)
    colors = [Accent8[c] for c in clusters]
    # TODO: selection tools
    source = ColumnDataSource(dict(
        x=reduced_data[0,:],
        y=reduced_data[1,:],
        color=colors,
        label=clusters))
    f1 = figure(tools='tap,wheel_zoom,pan,lasso_select')
    f1.circle(x='x', y='y',
        color='color', fill_alpha=0.7, size=5,
        legend='label', source=source)
    f1.xaxis.axis_label = 'dim1'
    f1.yaxis.axis_label = 'dim2'
    s2 = ColumnDataSource(data=dict(x=[], y=[]))
    source.callback = CustomJS(args=dict(s2=s2), code="""
            var inds = cb_obj.selected['1d'].indices;
            """)
    script, div = components(f1)

    # no need to have downloads for m and w since it can be assumed that the user already has m and w
    env = Environment(loader=PackageLoader('app', 'templates'))
    template = env.get_template('state_estimation_static.html')
    output = template.render(user_id='', has_result=True, visualization=div+'\n\n'+script)
    with open(output_file, 'w') as f:
        f.write(output)
