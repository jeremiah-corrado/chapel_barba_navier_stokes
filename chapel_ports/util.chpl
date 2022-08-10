use IO;
use IO.FormattedIO;
import FileSystem.mkdir;
import Path.splitPath;

proc write_array_to_file(path: string, a : [?D] real) {
    const (file_path, file_name) = splitPath(path);

    try {
        mkdir(file_path, parents = true);
    } catch {
        writeln("unable to make specified directory: ", file_path);
    }

    var f: file;
    try {
        f = open(path, iomode.cw);
        var write_channel = f.writer();
        write_array_to_channel(write_channel, a);
    } catch {
        writeln("Unable to write to file: ", file_name);
    }

    try {
        f.close();
    } catch {
        writeln("Unable to close file: ", file_name);
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

proc linspace_dist(min: real, max: real, n: int, D) {
    const step = (max - min) / (n - 1):real;
    var a : [D] real;
    forall i in 0..#n {
        a[i] = i * step;
    }
    return a;
}
