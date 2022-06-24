use IO;
use IO.FormattedIO;

proc write_array_to_file(path: string, a : [] real) {
    var f: file;
    try {
        f = open(path, iomode.cw);
        try {
            var write_channel = f.writer();
            try {
                // write_channel.writeln(a);
                for val in a {
                    write_channel.writef("%.8dr\n", val);
                }
            } catch {
                writeln("Write Failed!");
            }
        } catch {
            writeln("Unable to open channel to file: ", path);
        }
    } catch {
        writeln("Unable to open file: ", path);
    }

    try {
        f.close();
    } catch {
        writeln("Unable to close file: ", path);
    }
}

proc linspace(min: real, max: real, n: int) {
    const step = (max - min) / (n-1):real;
    var a : [{0..<n}] real;
    for i in 0..#n {
        a[i] = i * step;
    }
    return a;
}
