using Literate

"""
    notebook_filter(str)

Filter to correctly render inline math (two backquotes) in notebooks.
"""
function notebook_filter(str)
    re = r"(?<!`)``(?!`)"  # Two backquotes not preceded by nor followed by another
    replace(str, re => "\$")
end

"""
    notebook(input_file)

Create notebooks from `input_file`, applying needed filters and write output to
the root folder of our typical Jupyter session.
"""
function notebook(input_file)
    output_dir = "."
    Literate.notebook(input_file, output_dir; preprocess=notebook_filter)
end

nothing;

