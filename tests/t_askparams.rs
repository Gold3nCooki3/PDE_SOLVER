mod common;
extern crate assert_cli;

#[cfg(test)]
mod integration {
    use assert_cli;

    #[test]
    fn test_ask_params(){
        assert_cli::Assert::main_binary()
            .with_args(&["1","1","0","1","2","1"])
            .unwrap();
        assert_cli::Assert::main_binary()
            .with_args(&["1","1","0","1","1","1e-5"])
            .unwrap();
    }

    #[test]
     fn ask_params_input_values(){
        assert_cli::Assert::main_binary()
            .with_args(&["1","1","0","1","2","200000"])
            .fails()
            .unwrap();

        assert_cli::Assert::main_binary()
            .with_args(&["1","1","0","1","1","5e-21"])
            .fails()
            .unwrap();

        assert_cli::Assert::main_binary()
            .with_args(&["1","1","200000","1","2","1"])
            .fails()
            .unwrap();

        assert_cli::Assert::main_binary()
            .with_args(&["1","3","20","1","3","1"])
            .fails()
            .unwrap();
    }
    
    #[test]
    fn ask_params_input_wrong(){
        assert_cli::Assert::main_binary()
            .with_args(&["1","1","abc","1","2","1"])
            .fails()
            .unwrap();

        assert_cli::Assert::main_binary()
            .with_args(&["ab","1","0","1","2","d"])
            .fails()
            .unwrap();

        assert_cli::Assert::main_binary()
            .with_args(&["1","!","0","1","%",","])
            .fails()
            .unwrap();
    }

   
    #[test]
    fn ask_params_input_less(){
        assert_cli::Assert::main_binary()
            .with_args(&["1","1","0","1","2"])
            .fails()
            .unwrap();
        assert_cli::Assert::main_binary()
            .with_args(&["1","1","2"])
            .fails()
            .unwrap();
        assert_cli::Assert::main_binary()
            .fails()
            .unwrap();
    }

}