use IO;
use IO.FormattedIO;

proc write_array_to_file(path: string, a : [?D] real) {
    var f: file;
    try! {
        f = open(path, iomode.cw);
        var write_channel = f.writer();
        write_array_to_channel(write_channel, a);
    } catch {
        writeln("Unable to write to file: ", path);
    }

    try {
        f.close();
    } catch {
        writeln("Unable to close file: ", path);
    }
}

private proc write_array_to_channel(channel, a : [?D] real) where D.rank == 1 {
    for val in a {
        channel.writef("%.8dr\n", val);
    }
}

private proc write_array_to_channel(channel, a : [?D] real) where D.rank == 2 {
    const (rx, ry) = D.dims();
    for i in rx {
        for j in (ry.first)..(ry.last - 1) {
            channel.writef("%.8dr ", a[(i, j)]);
        }
        channel.writef("%.8dr\n", a[(i, ry.last)]);
    }
}

proc linspace(min: real, max: real, n: int) {
    const step = (max - min) / (n - 1):real;
    var a : [{0..<n}] real;
    for i in 0..#n {
        a[i] = i * step;
    }
    return a;
}
