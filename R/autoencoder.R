library(tidyverse)
library(keras)
source("R/utils_fit.R")

opar <- par(no.readonly = TRUE)

dir_step2 <- "intermediates/2502/250514_step2/"


# load ----
smooth_preds <- list.files(dir_step2,
                           pattern = "_preds\\.qs$") |>
  enframe(value = "filename",
          name = NULL) |>
  separate_wider_regex(filename,
                       patterns = c(
                         cell_type = "^.+",
                         "_preds\\.qs"
                       ),
                       cols_remove = FALSE) |>
  pmap(\(cell_type, filename){
    mat_preds <- qs::qread(file.path(dir_step2,
                                     filename))
    colnames(mat_preds) <- paste0(cell_type, "|", colnames(mat_preds))
    
    mat_preds
  }) |>
  do.call(cbind, args = _)


smooth_centered <- circ_perm_mat(smooth_preds)

len <- nrow(smooth_centered)




par(mfrow = c(1,2), mar = c(3, 2, 2, 1) + 0.1)
exple <- sample(ncol(smooth_preds), 20)
matplot(smooth_preds[,exple], type = "l")
matplot(smooth_centered[,exple], type = "l")
abline(v = len/2, lty = "dotted", col = "grey")
par(opar)


par(mfrow = c(1,2), mar = c(3, 2, 2, 1) + 0.1)
matplot(log1p(smooth_preds)[,exple], type = "l")
matplot(log1p(smooth_centered)[,exple], type = "l")
abline(v = len/2, lty = "dotted", col = "grey")
par(opar)



# autoencoder ----
genes <- colnames(smooth_centered)

x_train <- log1p(smooth_centered)[, sample(genes, .8*length(genes))] |> t()
x_test <- log1p(smooth_centered)[, setdiff(genes, rownames(x_train))] |> t()



stopifnot(all.equal(
  genes |> sort(),
  union(rownames(x_train), rownames(x_test)) |> sort()
))

dim(x_train);dim(x_test)

# matplot(x_train[sample(nrow(x_train), 5),] |> t(), type = "l")


tensorflow::tf$random$set_seed(123)
set.seed(123)

n_bottleneck <- 5

# Define the encoder
encoder <- keras_model_sequential() |>
  layer_dense(units = 256, input_shape = len, activation = 'relu') |>
  layer_dense(units = n_bottleneck, activation = NULL, name = "bottleneck") |>
  layer_activation_leaky_relu(alpha = .6)

# Define the decoder
decoder <- keras_model_sequential() |>
  layer_dense(units = 256, activation = 'relu', input_shape = n_bottleneck) |>
  layer_dense(units = len, activation = "linear")

# Connect them to create the autoencoder
autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))
autoencoder  |> compile(optimizer = optimizer_adam(learning_rate = 0.0005), loss = 'mean_squared_error')


fitted <- autoencoder |>
  fit(x_train,
      x_train,
      epochs = 10,
      batch_size = 32,
      validation_data = list(x_test, x_test),
      verbose = TRUE)


# extract data ----
intermediate_layer_model <- keras_model(inputs = autoencoder$input,
                                        outputs = get_layer(autoencoder, "bottleneck")$output)
intermediate_output <- predict(intermediate_layer_model, t(smooth_centered))


reconstructed <- predict(autoencoder, t(smooth_centered))


rownames(intermediate_output) <- rownames(reconstructed) <- colnames(smooth_centered)



# check result ----
interms <- scale(intermediate_output)

par(mar = c(2, 2, 2, 1) + 0.1)
layout(matrix(c(1,3,2,3), nrow = 2))
exple <- sample(rownames(x_train), 5)


matplot(t(x_train)[,exple], type = "l")
matplot(t(reconstructed)[,exple], type = "l")
matplot(t(interms)[,exple], type = "p")
par(opar)





#~ save ----

# qs::qsave(intermediate_output,
#           file.path(dir_step2, "250515_autoencoder.qs"))




