#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100
#define W 512
#define patch "save/save.txt"
typedef struct veci2 {
    int x, y;
}veci2_t, vectorInt2_t;
typedef struct vec2 {
    double x, y;
}vec2_t, vector2_t;
typedef struct vecf2 {
    float x, y;
}vecf2_t, vectorFloat2_t;
typedef struct barnes {
    vecf2_t cords; // координаты центра масс
    int mass; //масса узла
    struct barnes * a , * b , * c , * d ; // указатели на поддеревья 
    int rx2; // диаметр узла
}barnes_t;
typedef struct point{
    vecf2_t cords;
    vecf2_t speed;
    vecf2_t acceleration;
}point_t;
struct point *points_w;
barnes_t *barnes_w;
FILE *save;
float distancef(struct point*, barnes_t*, barnes_t*);
point_t* noise(int num_p, int w_size){
    //printf("%i",num_p);
    point_t *points = (point_t*)calloc(num_p, sizeof(point_t));
    for(int i = 0; i < num_p; i++){
        points[i].cords.x = rand()%w_size;
        points[i].cords.y = rand()%w_size;
        points[i].speed.x = 0;
        points[i].speed.y = 0;
        points[i].acceleration.x = 0;
        points[i].acceleration.y = 0;
    }
    return points;
}
point_t* init_points(point_t* (*op)(int, int), int num_p, int w_size){
   return op(num_p, w_size);
}
barnes_t* update_barnes(int count, int nx, int ny, int ex, int ey){
    barnes_t *barnes = (barnes_t*)malloc(sizeof(barnes_t));
    barnes->cords.x = 0;
    barnes->cords.y = 0;
    barnes->rx2 = ex - nx;
    int num = 0;
    for(int i = 0; i < N; i++){
        if(points_w[i].cords.x < ex && points_w[i].cords.y < ey &&
           points_w[i].cords.x > nx && points_w[i].cords.y > ny)
        {
            num++;
            barnes->cords.x += points_w[i].cords.x;
            barnes->cords.y += points_w[i].cords.y;
        }
    }
    if(num != 0){
        barnes->cords.x = barnes->cords.x / num;
        barnes->cords.y = barnes->cords.y / num;
    }
    else{
        barnes->cords.x = 0;
        barnes->cords.y = 0;
    }
    barnes->mass = num;
    if(num <= 1){
        //printf("%i\n",num);
        char str [50];
        snprintf(str, sizeof(str), "b %i %i %i %i\n", nx, ny, ex, ey);
        if(save){
        fputs(str, save);
        }
        return barnes;
    }
    else{
        //printf("%i\n",count);
        int nx2 = nx >> 1;
        int ny2 = ny >> 1;
        int ex2 = ex >> 1;
        int ey2 = ey >> 1;
        barnes->a = update_barnes(count + 1, nx, ny, nx2 + ex2, ny2 + ey2);
        barnes->b = update_barnes(count + 1, (ex - ey2) + ny2, ny, ex, ey2 + ny2);
        barnes->c = update_barnes(count + 1, nx2 + ex2, ny2 + ey2, ex, ey);
        barnes->d = update_barnes(count + 1, nx, ny2 + ey2, nx2 + ex2, ey); 
        return barnes;
    }
}
void dispoas_barnes(barnes_t* barnes_to_dispoas){
    if(barnes_to_dispoas->a != NULL){
        dispoas_barnes(barnes_to_dispoas->a);
        barnes_to_dispoas->a = NULL;
    }
    if(barnes_to_dispoas->b != NULL){
        dispoas_barnes(barnes_to_dispoas->b);
        barnes_to_dispoas->b = NULL;
    }
    if(barnes_to_dispoas->c != NULL){
        dispoas_barnes(barnes_to_dispoas->c);
        barnes_to_dispoas->c = NULL;
    }
    if(barnes_to_dispoas->d != NULL){
        dispoas_barnes(barnes_to_dispoas->d);
        barnes_to_dispoas->d = NULL;
    }
    free(barnes_to_dispoas);
    barnes_to_dispoas = NULL;
    return;
}
void update_points(void){
    for(int i = 0; i < N; i++){
    char str [50];
    snprintf(str, sizeof(str), "p %f %f\n", points_w[i].cords.x, points_w[i].cords.y);
    if(save){
        fputs(str, save);
    }
    barnes_t *barnes_locale;
    float distance = distancef(&points_w[i], barnes_locale, barnes_w);
    //printf("%f\n", distance);
    float aot = (barnes_locale->mass) / (distance * distance);
    points_w[i].acceleration.x = aot * ((barnes_locale->cords.x-points_w[i].cords.x)/distance);
    points_w[i].acceleration.y = aot * ((barnes_locale->cords.y-points_w[i].cords.y)/distance);
    //points_w[i].acceleration.x *= 0.0001f;
    //points_w[i].acceleration.y *= 0.0001f;
    points_w[i].speed.x += points_w[i].acceleration.x;
    points_w[i].speed.y += points_w[i].acceleration.y;
    points_w[i].cords.x += points_w[i].speed.x;
    points_w[i].cords.y += points_w[i].speed.y;
    }

}
float distancef(struct point*point, barnes_t *barnes_l,barnes_t *barnes_locale_w){
    //printf("%f %f %f %f\n", barnes_locale_w->cords.x, barnes_locale_w->cords.y, point->cords.x, point->cords.y);
    float x = (barnes_locale_w->cords.x - point->cords.x)*(barnes_locale_w->cords.x - point->cords.x);
    float y = (barnes_locale_w->cords.y - point->cords.y)*(barnes_locale_w->cords.y - point->cords.y);
    //printf("%f %f\n", x, y);
    float distance = sqrt(x + y);
    //printf("%f\n", distance);
    if(barnes_locale_w->rx2/distance < 0.5)
    {
        float a;
        float b;
        float c;
        float d;
        if(barnes_locale_w->a != NULL){
            //printf("a\n");
            a = distancef(point, barnes_l, barnes_locale_w->a);
            if(!isnan(a)){
                return a;
            }
        }
        if(barnes_locale_w->b != NULL){
            //printf("b\n");
            b = distancef(point, barnes_l, barnes_locale_w->b);
            if(!isnan(b)){
                return b;
            }
        }
        if(barnes_locale_w->c != NULL){
            //printf("c\n");
            c = distancef(point, barnes_l, barnes_locale_w->c);
            if(!isnan(c)){
                return c;
            }
        }
        if(barnes_locale_w->d != NULL){
            //printf("d\n");
            d = distancef(point, barnes_l, barnes_locale_w->d);
            if(!isnan(d)){
                return d;
            }
        }
    }
    else{
        //printf("else\n");
        barnes_l = barnes_locale_w;
        return distance;
    }
}
int main(int argc, char** argv)
{
    save = fopen(patch, "w"); 
    points_w = init_points(noise, N, W);
    int count = 0;
    while(count < 1000){
    barnes_w = update_barnes(1,0,0,W,W);
    //printf("x: %f y: %f\n",barnes_w->cords.x, barnes_w->cords.y);
    update_points();
    dispoas_barnes(barnes_w);
    char frame [20];
    snprintf(frame, sizeof(frame), "frame #%i\n", count);
    if(save){
        fputs(frame, save);
    }
    count++;
    }
    fclose(save); 
    //for(int i = 0; i < 1; i++){
    //printf("a:%i b:%i c:%i d:%i \n",barnes_w->a->mass,barnes_w->b->mass,barnes_w->c->mass,barnes_w->d->mass);
    //printf("a x: %i y: %i\nb x: %i y: %i\nc x: %i y: %i\nd x: %i y: %i\n",barnes_w->a->a->cords.x, barnes_w->a->a->cords.y,barnes_w->a->b->cords.x, barnes_w->a->b->cords.y,barnes_w->a->c->cords.x, barnes_w->a->c->cords.y,barnes_w->a->d->cords.x, barnes_w->a->d->cords.y);
    //printf("a x: %i y: %i\nb x: %i y: %i\nc x: %i y: %i\nd x: %i y: %i\n",barnes_w->c->a->cords.x, barnes_w->c->a->cords.y,barnes_w->c->b->cords.x, barnes_w->c->b->cords.y,barnes_w->c->c->cords.x, barnes_w->c->c->cords.y,barnes_w->c->d->cords.x, barnes_w->c->d->cords.y);
    //}
    return 0;
}
