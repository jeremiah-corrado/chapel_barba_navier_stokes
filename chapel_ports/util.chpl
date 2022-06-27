use IO;
use IO.FormattedIO;

proc write_array_to_file(path: string, a : [?D] real) {
    var f: file;
    try! {
        f = open(path, iomode.cw);
        var write_channel = f.writer();
        select D.rank {
            when 1 {
                for val in a {
                    write_channel.writef("%.8dr\n", val);
                }
            }
            when 2 {
                const (rx, ry) = D.dims();
                for i in rx {
                    for j in (ry.first)..(ry.last - 1) {
                        write_channel.writef("%.8dr ", a[(i, j)]);
                    }
                    write_channel.writef("%.8dr\n", a[(i, ry.last)]);
                }
            }
            otherwise {
                writeln("Only rank 1 and 2 arrays are supported!");
            }
        }
    } catch {
        writeln("Unable to write to file: ", path);
    }

    try {
        f.close();
    } catch {
        writeln("Unable to close file: ", path);
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
