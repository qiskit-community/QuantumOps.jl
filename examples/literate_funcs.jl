function notebook_filter(str)
    re = r"(?<!`)``(?!`)"  # Two backquotes not preceded by nor followed by another
    replace(str, re => "\$")
end

function notebook(input_file)
    Literate.notebook(input_file, "/home/lapeyre/code/github/Qiskit/oneenv/"; preprocess=notebook_filter)
end

