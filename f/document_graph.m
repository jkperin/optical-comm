function m = document_graph(h, filename)
    m = matlab2tikz(h);

    m.write(filename);
