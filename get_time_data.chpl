import Subprocess as Sub;
import List;
import Time;
import FileSystem as FS;
use IO;
use IO.FormattedIO;

config const n = 10,
             nl = 1,
             step = 11,
             dist = false;

const program_name =  "step_" + step:string + (if dist then "_dist" else "");

proc main(args: [] string) {
    writeln("Testing: ", program_name);

    // // compile: chpl --fast -o bin/{name} chapel_ports/{name}.chpl
    // writeln("Compiling...");
    // var comp = Sub.spawn(["chpl", "--fast", "-o", "bin/" + program_name, "./chapel_ports/" + program_name + ".chpl"]);
    // comp.wait();

    // run: ./bin/{name} -nl {nl} {args}
    writeln("Running...");
    var spawn_args_list = new List.list(["./bin/" + program_name, "-nl", nl:string]); spawn_args_list.append(args[1..]);
    const spawn_args_array = spawn_args_list.toArray();

    var walltimes : [0..<n] real;
    coforall test_id in 0..#n {
        var t = new Time.Timer();
        t.start();
        var prog = Sub.spawn(spawn_args_array);
        prog.wait();
        walltimes[test_id] = t.elapsed();
    }

    writeln("Walltimes (sec): ", walltimes);
    const mean = (+ reduce walltimes) / (n:real);
    writeln("Mean Walltime: ", mean);

    writeln("Exporting Data...");
    exportData("./perf_data_" + program_name + ".csv", walltimes, mean, spawn_args_array);
}

proc exportData(fileName: string, walltimes: [] real, mean: real, spawn_args: [] string) {
    var f, w, r;

    try! {
        if !FS.exists(fileName) then writeHeader(fileName);
        f = open(fileName, iomode.r);
        r = f.reader();

        var current_contents: string;
        r.readstring(current_contents);

        r.close();
        f.close();

        f = open(fileName, iomode.cw);
        w = f.writer();

        w.write(current_contents);
        w.write("".join(spawn_args));
        for wt in walltimes do w.writef(", %.5dr", wt);
        w.writef(", %.5dr\n", mean);

        w.close();
        f.close();
    }
    
}

proc writeHeader(fileName: string) {
    var f, w;

    try! { 
        f = open(fileName, iomode.cw);
        w = f.writer();

        w.write("Trial: ");
        for tid in 0..#n do w.write(", ", tid:string);
        w.writeln(", mean");

        w.close();
        f.close();
    }
}