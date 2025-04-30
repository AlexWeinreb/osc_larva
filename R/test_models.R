# interactive testing code


source("R/test_models.R")
to_test <- as.character(c(27,31:32))

models <- set_names(to_test, paste0("mod", to_test)) |>
  map(
    ~ get(paste0("ae_", .x))
  )

res <- map(1:10,
           \(.rep){
             
             imap(models,
                  \(.mod, .nm){
                    tibble(model = .nm,
                           val_loss = .mod())
                  }
             ) |>
               bind_rows() |>
               add_column(replicate = .rep)
             
           },
           .progress = TRUE) |>
  bind_rows()

res

ggplot(res) +
  ggbeeswarm::geom_quasirandom(aes(x = model, y = val_loss))

res |>
  ggplot() +
  # coord_cartesian(ylim = c(NA, .01)) +
  geom_boxplot(aes(x = model, y = val_loss))












ae_1 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_conv_1d(filters = 64, kernel_size = 3, activation = 'relu', input_shape = c(len, 1),
                  name = "inputConv") |>
    layer_max_pooling_1d(pool_size = 2,
                         name = "Pool") |>
    layer_flatten() |>
    layer_dense(units = n_bottleneck, activation = 'relu', name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 64, activation = 'relu', input_shape = n_bottleneck,
                name = "denseFromBottleneck") |>
    layer_reshape(target_shape = c(64, 1),
                  name = "reshape") |>
    layer_conv_1d_transpose(filters = 64, kernel_size = 3,
                            strides = 2, activation = 'relu',
                            padding = 'same',
                            name = "deconv") |>
    layer_conv_1d(filters = 1, kernel_size = 3,
                  activation = 'linear', padding = 'same',
                  name = "output")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 128,
        validation_data = list(x_test, x_test))
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}


ae_2 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_conv_1d(filters = 64, kernel_size = 5, activation = 'relu', input_shape = c(len, 1),
                  name = "inputConv") |>
    layer_max_pooling_1d(pool_size = 2,
                         name = "Pool") |>
    layer_flatten() |>
    layer_dense(units = n_bottleneck, activation = 'relu', name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 64, activation = 'relu', input_shape = n_bottleneck,
                name = "denseFromBottleneck") |>
    layer_reshape(target_shape = c(64, 1),
                  name = "reshape") |>
    layer_conv_1d_transpose(filters = 64, kernel_size = 5,
                            strides = 2, activation = 'relu',
                            padding = 'same',
                            name = "deconv") |>
    layer_conv_1d(filters = 1, kernel_size = 3,
                  activation = 'linear', padding = 'same',
                  name = "output")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 128,
        validation_data = list(x_test, x_test))
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}


ae_3 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_conv_1d(filters = 64, kernel_size = 7, activation = 'relu', input_shape = c(len, 1),
                  name = "inputConv") |>
    layer_max_pooling_1d(pool_size = 2,
                         name = "Pool") |>
    layer_flatten() |>
    layer_dense(units = n_bottleneck, activation = 'relu', name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 64, activation = 'relu', input_shape = n_bottleneck,
                name = "denseFromBottleneck") |>
    layer_reshape(target_shape = c(64, 1),
                  name = "reshape") |>
    layer_conv_1d_transpose(filters = 64, kernel_size = 7,
                            strides = 2, activation = 'relu',
                            padding = 'same',
                            name = "deconv") |>
    layer_conv_1d(filters = 1, kernel_size = 3,
                  activation = 'linear', padding = 'same',
                  name = "output")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 128,
        validation_data = list(x_test, x_test))
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}


ae_4 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_conv_1d(filters = 64,
                  kernel_size = 3,
                  activation = 'relu',
                  input_shape = c(len, 1),
                  name = "inputConv") |>
    layer_max_pooling_1d(pool_size = 2,
                         name = "Pool") |>
    layer_conv_1d(filters = 32,
                  kernel_size = 3,
                  activation = 'relu') |>
    layer_max_pooling_1d(pool_size = 2) |>
    layer_flatten() |>
    layer_dense(units = n_bottleneck,
                activation = 'relu',
                name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 32, activation = 'relu', input_shape = n_bottleneck,
                name = "denseFromBottleneck") |>
    layer_reshape(target_shape = c(32, 1),
                  name = "reshape") |>
    layer_conv_1d_transpose(filters = 32,
                            kernel_size = 3,
                            strides = 2,
                            activation = 'relu',
                            padding = 'same',
                            name = "deconv") |>
    layer_conv_1d_transpose(filters = 64,
                            kernel_size = 3,
                            strides = 2,
                            activation = 'relu',
                            padding = 'same') |>
    layer_conv_1d(filters = 1,
                  kernel_size = 3,
                  activation = 'linear',
                  padding = 'same',
                  name = "output")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 128,
        validation_data = list(x_test, x_test))
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}




ae_5 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_conv_1d(filters = 64,
                  kernel_size = 3,
                  activation = 'relu',
                  input_shape = c(len, 1),
                  name = "inputConv") |>
    layer_max_pooling_1d(pool_size = 2,
                         name = "Pool") |>
    layer_conv_1d(filters = 32,
                  kernel_size = 7,
                  activation = 'relu') |>
    layer_max_pooling_1d(pool_size = 2) |>
    layer_flatten() |>
    layer_dense(units = n_bottleneck,
                activation = 'relu',
                name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 32, activation = 'relu', input_shape = n_bottleneck,
                name = "denseFromBottleneck") |>
    layer_reshape(target_shape = c(32, 1),
                  name = "reshape") |>
    layer_conv_1d_transpose(filters = 32,
                            kernel_size = 7,
                            strides = 2,
                            activation = 'relu',
                            padding = 'same',
                            name = "deconv") |>
    layer_conv_1d_transpose(filters = 64,
                            kernel_size = 3,
                            strides = 2,
                            activation = 'relu',
                            padding = 'same') |>
    layer_conv_1d(filters = 1,
                  kernel_size = 3,
                  activation = 'linear',
                  padding = 'same',
                  name = "output")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 128,
        validation_data = list(x_test, x_test))
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}





ae_6 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 64,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = 'relu',
                name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 64, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 128,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}



ae_7 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = 'relu',
                name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 128,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}


ae_8 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = 'relu',
                name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 64,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}




ae_9 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 256,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = 'relu',
                name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 256, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 64,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}





ae_10 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = 'relu',
                name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 32,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}







ae_11 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 256,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = 'relu',
                name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 256, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 32,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}


ae_12 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = 'relu',
                name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}




ae_13 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 256,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = 'relu',
                name = "bottleneck")
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 256, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}




ae_14 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu()
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}




ae_15 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu()
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 32,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}



ae_16 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .2)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}


ae_17 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .1)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}





ae_18 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .4)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}



ae_19 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .4)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 32,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}






ae_20 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .5)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}



ae_21 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .6)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}



ae_22 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .7)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}

ae_23 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .8)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}

ae_24 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 128,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .9)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 128, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}


ae_25 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 256,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .6)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 256, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}




ae_26 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 64,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .6)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 64, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}




ae_27 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 512,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .6)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 512, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}


ae_28 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 256,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .6)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 256, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 32,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}

ae_29 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 1024,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .6)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 1024, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}



ae_30 <- function(){
  n_bottleneck <- 8
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 512,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .6)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 512, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 8,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}



ae_31 <- function(){
  n_bottleneck <- 6
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 512,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .6)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 512, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}



ae_32 <- function(){
  n_bottleneck <- 10
  
  # Define the encoder
  encoder <- keras_model_sequential() |>
    layer_dense(units = 512,
                input_shape = len,
                activation = 'relu') |>
    layer_dense(units = n_bottleneck,
                activation = NULL,
                name = "bottleneck") |>
    layer_activation_leaky_relu(alpha = .6)
  
  # Define the decoder
  decoder <- keras_model_sequential() |>
    layer_dense(units = 512, activation = 'relu', input_shape = n_bottleneck) |>
    layer_dense(units = len, activation = "linear")
  
  # Connect them to create the autoencoder
  autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
  autoencoder  |> compile(optimizer = 'adam', loss = 'mean_squared_error')
  
  
  fitted <- autoencoder |>
    fit(x_train,
        x_train,
        epochs = 20,
        batch_size = 16,
        validation_data = list(x_test, x_test),
        verbose = FALSE)
  
  
  fitted$metrics$val_loss[[ length(fitted$metrics$val_loss) ]]
  
}
