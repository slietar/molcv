fn main() {
    let args = std::env::args_os()
        .map(|arg| arg.into_string().unwrap())
        .collect::<Vec<_>>();

    molcv::cli(&args);
}
