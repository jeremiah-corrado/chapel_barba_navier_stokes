import Subprocess as Sub;
import List;
import Time;
import FileSystem as FS;
import Path;
use IO;
use IO.FormattedIO;

config const n = 10,
             nl = 1,
             step = 11,
             dist = false,
             recompile = false;

const program_name =  "step_" + step:string + (if dist then "_dist" else "");

proc main(args: [] string) {
    writeln("Testing: ", program_name);

    if recompile {
        // compile: chpl --fast -o bin/{name} chapel_ports/{name}.chpl
        writeln("Compiling...");
        var comp = Sub.spawn(["chpl", "--fast", "-o", "bin/" + program_name, "./chapel_ports/" + program_name + ".chpl"]);
        comp.wait();
    }

    const args_clone = args[1..];

    // run: ./bin/{name} -nl {nl} {args}
    writeln("Running Chapel...");
    var spawn_args_list = new List.list(["./bin/" + program_name, "-nl", nl:string]); spawn_args_list.append(args_clone);
    const spawn_args_array = spawn_args_list.toArray();
    const (chpl_walltimes, chpl_mean) = exec_tests(spawn_args_array);

    writeln("Exporting Chapel Data...");
    exportData("./perf_data_" + program_name + ".csv", chpl_walltimes, chpl_mean, spawn_args_array, "chpl");

    writeln("Running Python...");
    var py_spawn_args_list = new List.list([
        "srun", "--job-name=py_step_" + step:string, 
        "--quiet", "--nodes=1", "--ntasks=1", "--ntasks-per-node=1", "--cpus-per-task=96",
        "--exclusive", "--mem=0", "--kill-on-bad-exit", "--constraint=CL48,192GB", 
        "python3", Path.absPath("./python_scripts/step_" + step:string + ".py"):string
        ]);
    py_spawn_args_list.append(args[1..]);
    const py_spawn_args_array = py_spawn_args_list.toArray();
    const (py_walltimes, py_mean) = exec_tests(py_spawn_args_array);

    writeln("Exporting Python Data...");
    exportData("./perf_data_" + program_name + ".csv", py_walltimes, py_mean, spawn_args_array, "py");
}

proc exec_tests(spawn_args: [] string) {
    var walltimes : [0..<n] real;
    coforall test_id in 0..#n {
        var t = new Time.Timer();
        t.start();
        var prog = Sub.spawn(spawn_args);
        prog.wait();
        walltimes[test_id] = t.elapsed();
        assert(prog.exitCode == 0);
    }

    writeln("Walltimes (sec): ", walltimes);
    const mean = (+ reduce walltimes) / (n:real);
    writeln("Mean Walltime: ", mean);

    return (walltimes, mean);
}

proc exportData(fileName: string, walltimes: [] real, mean: real, spawn_args: [] string, impl: string) {
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
        w.write(impl, ", ");
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

        w.write("Impl, Trial: ");
        for tid in 0..#n do w.write(", ", tid:string);
        w.writeln(", mean");

        w.close();
        f.close();
    }
}
